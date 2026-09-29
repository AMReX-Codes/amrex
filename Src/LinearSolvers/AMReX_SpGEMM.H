#ifndef AMREX_SPGEMM_H_
#define AMREX_SPGEMM_H_
#include <AMReX_Config.H>

#include <AMReX_Algorithm.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_OpenMP.H>
#include <AMReX_SpMatrix.H>

#include <algorithm>
#include <limits>
#include <type_traits>
#include <utility>

namespace amrex::detail {

//! Row callback of SpGEMM that leaves the rows as they are.
struct NoRowPost {
    template <typename I, typename T>
    Long operator() (Long, I*, T*, Long n) const { return n; }
};

//! CSR with `nrows` empty rows.
template <typename T, template<typename> class V, typename I = Long>
CSR<T,V,I> spgemm_empty (Long nrows)
{
    CSR<T,V,I> C;
    C.row_offset.resize(nrows+1);
    auto* p = C.row_offset.data();
    ParallelForOMP(nrows+1, [=] AMREX_GPU_DEVICE (Long i) { p[i] = 0; });
    Gpu::streamSynchronize();
    return C;
}

#if !defined(AMREX_USE_GPU)

// Sorts the n entries of a row by column; `tmp` is scratch for long rows.
template <typename T, typename I>
void sort_row_cpu (I* AMREX_RESTRICT col, T* AMREX_RESTRICT val, Long n,
                   Vector<std::pair<I,T>>& tmp)
{
    if (n <= 32) {
        for (Long k = 1; k < n; ++k) {
            I const c = col[k];
            T const v = val[k];
            Long m = k;
            for (; m > 0 && col[m-1] > c; --m) {
                col[m] = col[m-1];
                val[m] = val[m-1];
            }
            col[m] = c;
            val[m] = v;
        }
    } else {
        tmp.resize(n);
        for (Long k = 0; k < n; ++k) { tmp[k] = {col[k], val[k]}; }
        std::sort(tmp.begin(), tmp.end(),
                  [] (auto const& x, auto const& y) { return x.first < y.first; });
        for (Long k = 0; k < n; ++k) { col[k] = tmp[k].first; val[k] = tmp[k].second; }
    }
}

/**
 * \brief Builds a CSR on the host one row at a time.
 *
 * `make_row(t)` returns a callable `row(i, col, val)` for block t that
 * computes row i in its own scratch storage, points `col` and `val` to it
 * and returns the number of entries. The rows are split into one block per
 * thread with about the same total `work(i)`. Rows are appended to buffers
 * and copied into the result once all counts are known, so no symbolic
 * pass is needed.
 */
template <typename T, template<typename> class V, typename I, typename W, typename F>
CSR<T,V,I> csr_from_rows_cpu (Long nrows, W const& work, F const& make_row)
{
    CSR<T,V,I> C;
    C.row_offset.resize(nrows+1);
    I* AMREX_RESTRICT crow = C.row_offset.data();
    AMREX_ASSUME(crow != nullptr); // gcc -Wnull-dereference false positive
    crow[0] = 0;

    int const nblocks = int(std::max(Long(1), std::min(Long(OpenMP::get_max_threads()), nrows)));
    Vector<Long> rbegin(nblocks+1, nrows);
    Long* AMREX_RESTRICT prb = rbegin.data();
    AMREX_ASSUME(prb != nullptr); // gcc -Wnull-dereference false positive
    prb[0] = 0;
    if (nblocks > 1) {
        Vector<Long> wsum(nrows+1);
        Long* AMREX_RESTRICT pw = wsum.data();
        AMREX_ASSUME(pw != nullptr);
        pw[0] = 0;
        for (Long i = 0; i < nrows; ++i) { pw[i+1] = pw[i] + work(i); }
        for (int t = 1; t < nblocks; ++t) {
            prb[t] = std::lower_bound(pw, pw+nrows+1, pw[nrows]/nblocks*t) - pw;
        }
    }

    if (nblocks == 1) {
        // Rows go straight into the result. The capacity is a guess from
        // `work`; pages that are never written are never touched.
        Long cap = 0;
        for (Long i = 0; i < nrows; ++i) { cap += work(i); }
        C.col_index.resize(cap);
        C.mat.resize(cap);
        auto row = make_row(0);
        Long p = 0;
        for (Long i = 0; i < nrows; ++i) {
            I const* rc = nullptr;
            T const* rv = nullptr;
            Long const n = row(i, rc, rv);
            if (p + n > Long(C.col_index.size())) {
                Long const newcap = std::max(p + n, Long(C.col_index.size())*3/2);
                C.col_index.resize(p); // only the filled part is copied
                C.mat.resize(p);
                C.col_index.resize(newcap);
                C.mat.resize(newcap);
            }
            std::copy(rc, rc+n, C.col_index.data()+p);
            std::copy(rv, rv+n, C.mat.data()+p);
            p += n;
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(p <= Long(std::numeric_limits<I>::max()),
                                             "SpGEMM: too many nonzeros for the index type");
            crow[i+1] = I(p);
        }
        C.col_index.resize(p);
        C.mat.resize(p);
        C.nnz = p;
        return C;
    }

    struct Buffer {
        V<I> col;
        V<T> val;
        Long n = 0;
    };
    Vector<Vector<Buffer>> buffers(nblocks);
    constexpr Long buffer_size = Long(1) << 20;

#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(static,1)
#endif
    for (int t = 0; t < nblocks; ++t) {
        auto row = make_row(t);
        auto& bufs = buffers[t];
        for (Long i = rbegin[t]; i < rbegin[t+1]; ++i) {
            I const* rc = nullptr;
            T const* rv = nullptr;
            Long const n = row(i, rc, rv);
            if (bufs.empty() || bufs.back().n + n > Long(bufs.back().col.size())) {
                bufs.emplace_back();
                bufs.back().col.resize(std::max(buffer_size, n));
                bufs.back().val.resize(std::max(buffer_size, n));
            }
            auto& buf = bufs.back();
            std::copy(rc, rc+n, buf.col.data()+buf.n);
            std::copy(rv, rv+n, buf.val.data()+buf.n);
            buf.n += n;
            crow[i+1] = I(n);
        }
    }

    // Exclusive scan in Long: the result may not fit the index type.
    Long total = 0;
    for (Long i = 0; i < nrows; ++i) {
        Long const cnt = crow[i+1];
        crow[i] = I(total);
        total += cnt;
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(total <= Long(std::numeric_limits<I>::max()),
                                     "SpGEMM: too many nonzeros for the index type");
    crow[nrows] = I(total);
    C.nnz = total;
    C.col_index.resize(C.nnz);
    C.mat.resize(C.nnz);
    I* AMREX_RESTRICT ccol = C.col_index.data();
    T* AMREX_RESTRICT cmat = C.mat.data();

#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(static,1)
#endif
    for (int t = 0; t < nblocks; ++t) {
        Long p = crow[rbegin[t]];
        for (auto& buf : buffers[t]) {
            std::copy(buf.col.data(), buf.col.data()+buf.n, ccol+p);
            std::copy(buf.val.data(), buf.val.data()+buf.n, cmat+p);
            p += buf.n;
            buf = Buffer{};
        }
    }

    return C;
}

// Gustavson's algorithm with a dense marker holding the position of each
// column in the current row. Output rows are sorted, then passed to
// `row_post(i, col, val, n)`, which may change them in place and returns
// the new length.
template <typename T, template<typename> class V, typename I, typename F = NoRowPost>
CSR<T,V,I> spgemm_local_cpu (Long nrows, Long ncols,
                             CsrView<T const,I> const& A, CsrView<T const,I> const& B,
                             F const& row_post = {})
{
    auto nprod = [&] (Long i) {
        Long n = 0;
        for (Long ap = A.row_offset[i]; ap < A.row_offset[i+1]; ++ap) {
            Long const k = A.col_index[ap];
            n += B.row_offset[k+1] - B.row_offset[k];
        }
        return n;
    };
    return csr_from_rows_cpu<T,V,I>(nrows, nprod, [&] (int) {
        return [&, marker = Vector<I>(ncols, I(-1)), col = Vector<I>(), val = Vector<T>(),
                tmp = Vector<std::pair<I,T>>()]
            (Long i, I const*& rc, T const*& rv) mutable -> Long
        {
            Long const maxlen = std::min(nprod(i), ncols);
            if (Long(col.size()) < maxlen) {
                col.resize(maxlen);
                val.resize(maxlen);
            }
            Long n = 0;
            for (Long ap = A.row_offset[i]; ap < A.row_offset[i+1]; ++ap) {
                Long const k = A.col_index[ap];
                T const a = A.mat[ap];
                for (Long bp = B.row_offset[k]; bp < B.row_offset[k+1]; ++bp) {
                    I const j = B.col_index[bp];
                    Long const m = marker[j]; // may be left over from an earlier row
                    if (m >= 0 && m < n && col[m] == j) {
                        val[m] += a * B.mat[bp];
                    } else {
                        marker[j] = I(n);
                        col[n] = j;
                        val[n] = a * B.mat[bp];
                        ++n;
                    }
                }
            }
            sort_row_cpu(col.data(), val.data(), n, tmp);
            n = row_post(i, col.data(), val.data(), n);
            rc = col.data();
            rv = val.data();
            return n;
        };
    });
}

#ifdef AMREX_USE_MPI
/**
 * \brief Local part of a distributed product on the host, without first
 * copying A and B into one block each.
 *
 * Row i of A is A0's row (column k) followed by A1's row (column nb + j).
 * Row k < nb of B is B0's row followed by B1's row with columns mapped by
 * `b1map`; row k >= nb is row k - nb of E. The sums are the same as those
 * of spgemm_local_cpu on the concatenated blocks. `row_post` is as in
 * spgemm_local_cpu.
 */
template <typename T, template<typename> class V, typename I, typename F>
CSR<T,V,I> spgemm_split_cpu (Long nrows, Long ncols, Long nb,
                             CsrView<T const,I> const& A0, CsrView<T const,I> const& A1,
                             CsrView<T const,I> const& B0, CsrView<T const,I> const& B1,
                             I const* b1map, CsrView<T const,I> const& E, F const& row_post)
{
    auto blen = [&] (Long k) -> Long {
        if (k < nb) {
            return (B0.row_offset[k+1] - B0.row_offset[k])
                + ((B1.nnz > 0) ? (B1.row_offset[k+1] - B1.row_offset[k]) : 0);
        } else {
            return E.row_offset[k-nb+1] - E.row_offset[k-nb];
        }
    };
    auto nprod = [&] (Long i) {
        Long n = 0;
        for (Long ap = A0.row_offset[i]; ap < A0.row_offset[i+1]; ++ap) { n += blen(A0.col_index[ap]); }
        if (A1.nnz > 0) {
            for (Long ap = A1.row_offset[i]; ap < A1.row_offset[i+1]; ++ap) { n += blen(nb + A1.col_index[ap]); }
        }
        return n;
    };
    return csr_from_rows_cpu<T,V,I>(nrows, nprod, [&] (int) {
        return [&, marker = Vector<I>(ncols, I(-1)), col = Vector<I>(), val = Vector<T>(),
                tmp = Vector<std::pair<I,T>>()]
            (Long i, I const*& rc, T const*& rv) mutable -> Long
        {
            Long const maxlen = std::min(nprod(i), ncols);
            if (Long(col.size()) < maxlen) {
                col.resize(maxlen);
                val.resize(maxlen);
            }
            Long n = 0;
            auto add = [&] (I j, T x) {
                Long const m = marker[j]; // may be left over from an earlier row
                if (m >= 0 && m < n && col[m] == j) {
                    val[m] += x;
                } else {
                    marker[j] = I(n);
                    col[n] = j;
                    val[n] = x;
                    ++n;
                }
            };
            auto brow = [&] (Long k, T a) {
                if (k < nb) {
                    for (Long bp = B0.row_offset[k]; bp < B0.row_offset[k+1]; ++bp) {
                        add(B0.col_index[bp], a * B0.mat[bp]);
                    }
                    if (B1.nnz > 0) {
                        for (Long bp = B1.row_offset[k]; bp < B1.row_offset[k+1]; ++bp) {
                            add(b1map[B1.col_index[bp]], a * B1.mat[bp]);
                        }
                    }
                } else {
                    for (Long bp = E.row_offset[k-nb]; bp < E.row_offset[k-nb+1]; ++bp) {
                        add(E.col_index[bp], a * E.mat[bp]);
                    }
                }
            };
            for (Long ap = A0.row_offset[i]; ap < A0.row_offset[i+1]; ++ap) {
                brow(A0.col_index[ap], A0.mat[ap]);
            }
            if (A1.nnz > 0) {
                for (Long ap = A1.row_offset[i]; ap < A1.row_offset[i+1]; ++ap) {
                    brow(nb + A1.col_index[ap], A1.mat[ap]);
                }
            }
            sort_row_cpu(col.data(), val.data(), n, tmp);
            n = row_post(i, col.data(), val.data(), n);
            rc = col.data();
            rv = val.data();
            return n;
        };
    });
}
#endif

/**
 * \brief Rows of A P on the host.
 *
 * Row i of A is A0's row (local column k) followed by A1's row (remote
 * index k). Row k of P is P0's row followed by P1's row with columns mapped
 * by `p1map`, and remote row k is row k of PE. `row(i, ...)` appends row i
 * to col/val at position p, with the same sums as spgemm_local_cpu, but
 * leaves it unsorted.
 */
template <typename T, typename I>
struct APRowsCpu
{
    CsrView<T const,I> A0, A1, P0, P1, PE;
    I const* p1map = nullptr;

    [[nodiscard]] Long plen (Long k) const {
        return (P0.row_offset[k+1] - P0.row_offset[k])
            + ((P1.nnz > 0) ? (P1.row_offset[k+1] - P1.row_offset[k]) : 0);
    }

    //! Upper bound on the length of row i.
    [[nodiscard]] Long maxlen (Long i) const {
        Long n = 0;
        for (Long ap = A0.row_offset[i]; ap < A0.row_offset[i+1]; ++ap) { n += plen(A0.col_index[ap]); }
        if (A1.nnz > 0) {
            for (Long ap = A1.row_offset[i]; ap < A1.row_offset[i+1]; ++ap) {
                Long const k = A1.col_index[ap];
                n += PE.row_offset[k+1] - PE.row_offset[k];
            }
        }
        return n;
    }

    //! `marker[j]` holds the position of column j; positions before p are
    //! from earlier rows. col/val must have room for maxlen(i) more entries.
    Long row (Long i, I* marker, I* col, T* val, Long p) const {
        Long const rs = p;
        auto add = [&] (I j, T x) {
            Long const m = marker[j]; // may be left over from elsewhere
            if (m >= rs && m < p && col[m] == j) {
                val[m] += x;
            } else {
                marker[j] = I(p);
                col[p] = j;
                val[p] = x;
                ++p;
            }
        };
        for (Long ap = A0.row_offset[i]; ap < A0.row_offset[i+1]; ++ap) {
            Long const k = A0.col_index[ap];
            T const a = A0.mat[ap];
            for (Long pp = P0.row_offset[k]; pp < P0.row_offset[k+1]; ++pp) {
                add(P0.col_index[pp], a * P0.mat[pp]);
            }
            if (P1.nnz > 0) {
                for (Long pp = P1.row_offset[k]; pp < P1.row_offset[k+1]; ++pp) {
                    add(p1map[P1.col_index[pp]], a * P1.mat[pp]);
                }
            }
        }
        if (A1.nnz > 0) {
            for (Long ap = A1.row_offset[i]; ap < A1.row_offset[i+1]; ++ap) {
                Long const k = A1.col_index[ap];
                T const a = A1.mat[ap];
                for (Long pp = PE.row_offset[k]; pp < PE.row_offset[k+1]; ++pp) {
                    add(PE.col_index[pp], a * PE.mat[pp]);
                }
            }
        }
        return p;
    }
};

/**
 * \brief Local part of the Galerkin product R (A P) on the host.
 *
 * Row c of R is R0's row (local fine rows) followed by R1's row (remote
 * index j, whose row of A P is row j of APE). The local rows of A P
 * (`ap`) are computed in chunks of fine rows when a row of R first needs
 * them, and a chunk is reused once no later row of R does, so A P is never
 * stored in full. The sums are done in the same order as
 * SpGEMM(R, SpGEMM(A, P)). Output rows are sorted.
 */
template <typename T, template<typename> class V, typename I>
CSR<T,V,I> rap_local_cpu (Long nf, Long nc, Long ncols,
                          CsrView<T const,I> const& R0, CsrView<T const,I> const& R1,
                          CsrView<T const,I> const& APE, APRowsCpu<T,I> const& ap)
{
    constexpr Long chunk_rows = 1024;
    Long const nchunks = (nf + chunk_rows - 1) / chunk_rows;

    // About the length of row c: also the capacity guess of csr_from_rows_cpu.
    auto work = [&] (Long c) {
        Long n = 0;
        for (Long rp = R0.row_offset[c]; rp < R0.row_offset[c+1]; ++rp) {
            Long const i = R0.col_index[rp];
            n += ap.A0.row_offset[i+1] - ap.A0.row_offset[i];
            if (ap.A1.nnz > 0) { n += ap.A1.row_offset[i+1] - ap.A1.row_offset[i]; }
        }
        if (R1.nnz > 0) {
            for (Long rp = R1.row_offset[c]; rp < R1.row_offset[c+1]; ++rp) {
                Long const j = R1.col_index[rp];
                n += APE.row_offset[j+1] - APE.row_offset[j];
            }
        }
        return n;
    };

    // Last row of R that needs a chunk: the largest local column of P in it.
    Vector<Long> chunk_last(nchunks, -1);
    for (Long q = 0; q < nchunks; ++q) {
        for (Long i = q*chunk_rows; i < std::min(nf, (q+1)*chunk_rows); ++i) {
            if (ap.P0.row_offset[i+1] > ap.P0.row_offset[i]) {
                chunk_last[q] = std::max(chunk_last[q],
                                         Long(ap.P0.col_index[ap.P0.row_offset[i+1]-1]));
            }
        }
    }

    struct Chunk {
        Vector<Long> off;
        Vector<I> col;
        Vector<T> val;
    };

    return csr_from_rows_cpu<T,V,I>(nc, work, [&] (int) {
        return [&, cmarker = Vector<I>(ncols, I(-1)), marker = Vector<I>(ncols, I(-1)),
                slot = Vector<int>(nchunks, -1), chunks = Vector<Chunk>(),
                live = Vector<Long>(), free_slots = Vector<int>(),
                next_release = std::numeric_limits<Long>::max(),
                col = Vector<I>(), val = Vector<T>(), tmp = Vector<std::pair<I,T>>()]
            (Long c, I const*& rc, T const*& rv) mutable -> Long
        {
            auto compute_chunk = [&] (Long q) {
                int s;
                if (free_slots.empty()) {
                    s = int(chunks.size());
                    chunks.emplace_back();
                } else {
                    s = free_slots.back();
                    free_slots.pop_back();
                }
                slot[q] = s;
                live.push_back(q);
                next_release = std::min(next_release, chunk_last[q]);
                auto& ch = chunks[s];
                Long const i0 = q*chunk_rows;
                Long const i1 = std::min(nf, i0+chunk_rows);
                ch.off.resize(i1-i0+1);
                ch.off[0] = 0;
                Long p = 0;
                for (Long i = i0; i < i1; ++i) {
                    Long const maxlen = ap.maxlen(i);
                    if (Long(ch.col.size()) < p + maxlen) {
                        ch.col.resize(std::max(p + maxlen, Long(ch.col.size())*2));
                        ch.val.resize(ch.col.size());
                    }
                    p = ap.row(i, cmarker.data(), ch.col.data(), ch.val.data(), p);
                    ch.off[i-i0+1] = p;
                }
            };

            // Chunks no later row needs are recycled.
            if (c > next_release) {
                next_release = std::numeric_limits<Long>::max();
                for (Long l = 0; l < Long(live.size()); ) {
                    Long const q = live[l];
                    if (chunk_last[q] < c) {
                        free_slots.push_back(slot[q]);
                        slot[q] = -1;
                        live[l] = live.back();
                        live.pop_back();
                    } else {
                        next_release = std::min(next_release, chunk_last[q]);
                        ++l;
                    }
                }
            }

            Long maxlen = 0;
            for (Long rp = R0.row_offset[c]; rp < R0.row_offset[c+1]; ++rp) {
                Long const i = R0.col_index[rp];
                Long const q = i / chunk_rows;
                if (slot[q] < 0) { compute_chunk(q); }
                auto const& ch = chunks[slot[q]];
                maxlen += ch.off[i-q*chunk_rows+1] - ch.off[i-q*chunk_rows];
            }
            if (R1.nnz > 0) {
                for (Long rp = R1.row_offset[c]; rp < R1.row_offset[c+1]; ++rp) {
                    Long const j = R1.col_index[rp];
                    maxlen += APE.row_offset[j+1] - APE.row_offset[j];
                }
            }
            maxlen = std::min(maxlen, ncols);
            if (Long(col.size()) < maxlen) {
                col.resize(maxlen);
                val.resize(maxlen);
            }

            Long n = 0;
            auto add = [&] (I j, T x) {
                Long const m = marker[j]; // may be left over from an earlier row
                if (m >= 0 && m < n && col[m] == j) {
                    val[m] += x;
                } else {
                    marker[j] = I(n);
                    col[n] = j;
                    val[n] = x;
                    ++n;
                }
            };
            for (Long rp = R0.row_offset[c]; rp < R0.row_offset[c+1]; ++rp) {
                Long const i = R0.col_index[rp];
                T const r = R0.mat[rp];
                Long const q = i / chunk_rows;
                auto const& ch = chunks[slot[q]];
                for (Long x = ch.off[i-q*chunk_rows]; x < ch.off[i-q*chunk_rows+1]; ++x) {
                    add(ch.col[x], r * ch.val[x]);
                }
            }
            if (R1.nnz > 0) {
                for (Long rp = R1.row_offset[c]; rp < R1.row_offset[c+1]; ++rp) {
                    Long const j = R1.col_index[rp];
                    T const r = R1.mat[rp];
                    for (Long x = APE.row_offset[j]; x < APE.row_offset[j+1]; ++x) {
                        add(APE.col_index[x], r * APE.mat[x]);
                    }
                }
            }
            sort_row_cpu(col.data(), val.data(), n, tmp);
            rc = col.data();
            rv = val.data();
            return n;
        };
    });
}

#elif defined(AMREX_USE_CUDA)

template <typename T, template<typename> class V, typename I>
CSR<T,V,I> spgemm_local_cusparse (Long nrows, Long ncols,
                                  CsrView<T const,I> const& A, CsrView<T const,I> const& B)
{
    static_assert(std::is_same_v<I,int>, "spgemm_local_cusparse: 32-bit indices only");

    cusparseHandle_t handle;
    AMREX_CUSPARSE_SAFE_CALL(cusparseCreate(&handle));
    AMREX_CUSPARSE_SAFE_CALL(cusparseSetStream(handle, Gpu::gpuStream()));

    cudaDataType data_type;
    if constexpr (std::is_same_v<T,float>) {
        data_type = CUDA_R_32F;
    } else if constexpr (std::is_same_v<T,double>) {
        data_type = CUDA_R_64F;
    } else if constexpr (std::is_same_v<T,GpuComplex<float>>) {
        data_type = CUDA_C_32F;
    } else if constexpr (std::is_same_v<T,GpuComplex<double>>) {
        data_type = CUDA_C_64F;
    } else {
        amrex::Abort("SpGEMM: unsupported data type");
    }

    CSR<T,V,I> C;
    C.row_offset.resize(nrows+1);

    AMREX_ALWAYS_ASSERT(nrows < Long(std::numeric_limits<int>::max()) &&
                        ncols < Long(std::numeric_limits<int>::max()));
    constexpr cusparseIndexType_t index_type = CUSPARSE_INDEX_32I;
    void* rowA = (void*)A.row_offset;
    void* colA = (void*)A.col_index;
    void* rowB = (void*)B.row_offset;
    void* colB = (void*)B.col_index;
    void* rowC = (void*)C.row_offset.data();

    cusparseSpMatDescr_t mat_A, mat_B, mat_C;
    AMREX_CUSPARSE_SAFE_CALL
        (cusparseCreateCsr(&mat_A, nrows, B.nrows, A.nnz, rowA, colA, (void*)A.mat,
                           index_type, index_type, CUSPARSE_INDEX_BASE_ZERO, data_type));
    AMREX_CUSPARSE_SAFE_CALL
        (cusparseCreateCsr(&mat_B, B.nrows, ncols, B.nnz, rowB, colB, (void*)B.mat,
                           index_type, index_type, CUSPARSE_INDEX_BASE_ZERO, data_type));
    AMREX_CUSPARSE_SAFE_CALL
        (cusparseCreateCsr(&mat_C, nrows, ncols, 0, rowC, nullptr, nullptr,
                           index_type, index_type, CUSPARSE_INDEX_BASE_ZERO, data_type));

    cusparseSpGEMMDescr_t spgemm_descr;
    AMREX_CUSPARSE_SAFE_CALL(cusparseSpGEMM_createDescr(&spgemm_descr));

    T alpha = T(1);
    T beta  = T(0);
    cusparseOperation_t op = CUSPARSE_OPERATION_NON_TRANSPOSE;

    // TODO: fall back to the chunked algorithm (CUSPARSE_SPGEMM_ALG3) when
    // the work buffer does not fit in memory.
    std::size_t buffer_size1 = 0;
    AMREX_CUSPARSE_SAFE_CALL
        (cusparseSpGEMM_workEstimation(handle, op, op, &alpha, mat_A, mat_B, &beta, mat_C,
                                       data_type, CUSPARSE_SPGEMM_DEFAULT, spgemm_descr,
                                       &buffer_size1, nullptr));
    auto* buffer1 = (void*)The_Arena()->alloc(buffer_size1);
    AMREX_CUSPARSE_SAFE_CALL
        (cusparseSpGEMM_workEstimation(handle, op, op, &alpha, mat_A, mat_B, &beta, mat_C,
                                       data_type, CUSPARSE_SPGEMM_DEFAULT, spgemm_descr,
                                       &buffer_size1, buffer1));

    std::size_t buffer_size2 = 0;
    AMREX_CUSPARSE_SAFE_CALL
        (cusparseSpGEMM_compute(handle, op, op, &alpha, mat_A, mat_B, &beta, mat_C,
                                data_type, CUSPARSE_SPGEMM_DEFAULT, spgemm_descr,
                                &buffer_size2, nullptr));
    auto* buffer2 = (void*)The_Arena()->alloc(buffer_size2);
    AMREX_CUSPARSE_SAFE_CALL
        (cusparseSpGEMM_compute(handle, op, op, &alpha, mat_A, mat_B, &beta, mat_C,
                                data_type, CUSPARSE_SPGEMM_DEFAULT, spgemm_descr,
                                &buffer_size2, buffer2));

    std::int64_t c_nrows, c_ncols, c_nnz;
    AMREX_CUSPARSE_SAFE_CALL(cusparseSpMatGetSize(mat_C, &c_nrows, &c_ncols, &c_nnz));
    AMREX_ALWAYS_ASSERT(c_nrows == nrows && c_ncols == ncols);
    AMREX_ALWAYS_ASSERT(c_nnz < std::int64_t(std::numeric_limits<int>::max()));

    C.mat.resize(c_nnz);
    C.col_index.resize(c_nnz);
    C.nnz = c_nnz;
    void* colC = (void*)C.col_index.data();
    AMREX_CUSPARSE_SAFE_CALL(cusparseCsrSetPointers(mat_C, rowC, colC, (void*)C.mat.data()));

    // cuSPARSE guarantees sorted column indices in C.
    AMREX_CUSPARSE_SAFE_CALL
        (cusparseSpGEMM_copy(handle, op, op, &alpha, mat_A, mat_B, &beta, mat_C,
                             data_type, CUSPARSE_SPGEMM_DEFAULT, spgemm_descr));

    Gpu::streamSynchronize();
    AMREX_CUSPARSE_SAFE_CALL(cusparseSpGEMM_destroyDescr(spgemm_descr));
    AMREX_CUSPARSE_SAFE_CALL(cusparseDestroySpMat(mat_A));
    AMREX_CUSPARSE_SAFE_CALL(cusparseDestroySpMat(mat_B));
    AMREX_CUSPARSE_SAFE_CALL(cusparseDestroySpMat(mat_C));
    AMREX_CUSPARSE_SAFE_CALL(cusparseDestroy(handle));
    The_Arena()->free(buffer1);
    The_Arena()->free(buffer2);

    return C;
}

#elif defined(AMREX_USE_HIP)

template <typename T, template<typename> class V, typename I>
CSR<T,V,I> spgemm_local_rocsparse (Long nrows, Long ncols,
                                   CsrView<T const,I> const& A, CsrView<T const,I> const& B)
{
    static_assert(std::is_same_v<I,int>, "spgemm_local_rocsparse: 32-bit indices only");

    rocsparse_handle handle;
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_create_handle(&handle));
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_set_stream(handle, Gpu::gpuStream()));

    AMREX_ALWAYS_ASSERT(nrows < Long(std::numeric_limits<int>::max()) &&
                        ncols < Long(std::numeric_limits<int>::max()));

    rocsparse_datatype data_type;
    if constexpr (std::is_same_v<T,float>) {
        data_type = rocsparse_datatype_f32_r;
    } else if constexpr (std::is_same_v<T,double>) {
        data_type = rocsparse_datatype_f64_r;
    } else if constexpr (std::is_same_v<T,GpuComplex<float>>) {
        data_type = rocsparse_datatype_f32_c;
    } else if constexpr (std::is_same_v<T,GpuComplex<double>>) {
        data_type = rocsparse_datatype_f64_c;
    } else {
        amrex::Abort("SpGEMM: unsupported data type");
    }

    constexpr rocsparse_indextype index_type = rocsparse_indextype_i32;
    constexpr rocsparse_index_base index_base = rocsparse_index_base_zero;

    CSR<T,V,I> C;
    C.row_offset.resize(nrows+1);

    rocsparse_spmat_descr mat_A, mat_B, mat_C, mat_D;
    AMREX_ROCSPARSE_SAFE_CALL
        (rocsparse_create_csr_descr(&mat_A, nrows, B.nrows, A.nnz,
                                    (void*)A.row_offset, (void*)A.col_index, (void*)A.mat,
                                    index_type, index_type, index_base, data_type));
    AMREX_ROCSPARSE_SAFE_CALL
        (rocsparse_create_csr_descr(&mat_B, B.nrows, ncols, B.nnz,
                                    (void*)B.row_offset, (void*)B.col_index, (void*)B.mat,
                                    index_type, index_type, index_base, data_type));
    AMREX_ROCSPARSE_SAFE_CALL
        (rocsparse_create_csr_descr(&mat_C, nrows, ncols, 0,
                                    (void*)C.row_offset.data(), nullptr, nullptr,
                                    index_type, index_type, index_base, data_type));
    // D is unused because beta is zero, but a valid descriptor is required.
    AMREX_ROCSPARSE_SAFE_CALL
        (rocsparse_create_csr_descr(&mat_D, 0, 0, 0, nullptr, nullptr, nullptr,
                                    index_type, index_type, index_base, data_type));

    T alpha = T(1);
    T beta  = T(0);
    auto const op = rocsparse_operation_none;
    auto const alg = rocsparse_spgemm_alg_default;

    std::size_t buffer_size = 0;
    AMREX_ROCSPARSE_SAFE_CALL
        (rocsparse_spgemm(handle, op, op, &alpha, mat_A, mat_B, &beta, mat_D, mat_C,
                          data_type, alg, rocsparse_spgemm_stage_buffer_size,
                          &buffer_size, nullptr));
    auto* buffer = (void*)The_Arena()->alloc(buffer_size);

    // This stage fills C's row offsets.
    AMREX_ROCSPARSE_SAFE_CALL
        (rocsparse_spgemm(handle, op, op, &alpha, mat_A, mat_B, &beta, mat_D, mat_C,
                          data_type, alg, rocsparse_spgemm_stage_nnz,
                          &buffer_size, buffer));

    std::int64_t c_nrows, c_ncols, c_nnz;
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_spmat_get_size(mat_C, &c_nrows, &c_ncols, &c_nnz));
    AMREX_ALWAYS_ASSERT(c_nrows == nrows && c_ncols == ncols);
    AMREX_ALWAYS_ASSERT(c_nnz < std::int64_t(std::numeric_limits<int>::max()));

    C.mat.resize(c_nnz);
    C.col_index.resize(c_nnz);
    C.nnz = c_nnz;
    AMREX_ROCSPARSE_SAFE_CALL
        (rocsparse_csr_set_pointers(mat_C, (void*)C.row_offset.data(),
                                    (void*)C.col_index.data(), (void*)C.mat.data()));

    AMREX_ROCSPARSE_SAFE_CALL
        (rocsparse_spgemm(handle, op, op, &alpha, mat_A, mat_B, &beta, mat_D, mat_C,
                          data_type, alg, rocsparse_spgemm_stage_compute,
                          &buffer_size, buffer));

    Gpu::streamSynchronize();
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_destroy_spmat_descr(mat_A));
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_destroy_spmat_descr(mat_B));
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_destroy_spmat_descr(mat_C));
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_destroy_spmat_descr(mat_D));
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_destroy_handle(handle));
    The_Arena()->free(buffer);

    C.sort(); // rocSPARSE does not promise sorted rows for the generic API.

    return C;
}

#elif defined(AMREX_USE_SYCL)

template <typename T, template<typename> class V, typename I>
CSR<T,V,I> spgemm_local_onemkl (Long nrows, Long ncols,
                                CsrView<T const,I> const& A, CsrView<T const,I> const& B)
{
    auto& q = Gpu::Device::streamQueue();

    CSR<T,V,I> C;
    C.row_offset.resize(nrows+1);
    // oneMKL wants valid pointers for C before the nnz is known.
    V<I> dummy_col(1);
    V<T> dummy_mat(1);

    mkl::sparse::matrix_handle_t hA{}, hB{}, hC{};
    mkl::sparse::init_matrix_handle(&hA);
    mkl::sparse::init_matrix_handle(&hB);
    mkl::sparse::init_matrix_handle(&hC);

#if defined(INTEL_MKL_VERSION) && (INTEL_MKL_VERSION < 20250300)
    mkl::sparse::set_csr_data(q, hA, I(nrows), I(B.nrows), mkl::index_base::zero,
                              (I*)A.row_offset, (I*)A.col_index, (T*)A.mat);
    mkl::sparse::set_csr_data(q, hB, I(B.nrows), I(ncols), mkl::index_base::zero,
                              (I*)B.row_offset, (I*)B.col_index, (T*)B.mat);
    mkl::sparse::set_csr_data(q, hC, I(nrows), I(ncols), mkl::index_base::zero,
                              C.row_offset.data(), dummy_col.data(), dummy_mat.data());
#else
    mkl::sparse::set_csr_data(q, hA, I(nrows), I(B.nrows), I(A.nnz), mkl::index_base::zero,
                              (I*)A.row_offset, (I*)A.col_index, (T*)A.mat);
    mkl::sparse::set_csr_data(q, hB, I(B.nrows), I(ncols), I(B.nnz), mkl::index_base::zero,
                              (I*)B.row_offset, (I*)B.col_index, (T*)B.mat);
    mkl::sparse::set_csr_data(q, hC, I(nrows), I(ncols), I(0), mkl::index_base::zero,
                              C.row_offset.data(), dummy_col.data(), dummy_mat.data());
#endif

    mkl::sparse::matmat_descr_t descr = nullptr;
    mkl::sparse::init_matmat_descr(&descr);
    mkl::sparse::set_matmat_data(descr,
                                 mkl::sparse::matrix_view_descr::general,
                                 mkl::transpose::nontrans,
                                 mkl::sparse::matrix_view_descr::general,
                                 mkl::transpose::nontrans,
                                 mkl::sparse::matrix_view_descr::general);

    using req = mkl::sparse::matmat_request;
    auto* size_buf = (std::int64_t*)The_Pinned_Arena()->alloc(sizeof(std::int64_t));

    mkl::sparse::matmat(q, hA, hB, hC, req::get_work_estimation_buf_size, descr,
                        size_buf, nullptr, {}).wait();
    auto* buffer1 = (void*)The_Arena()->alloc(std::size_t(*size_buf));
    mkl::sparse::matmat(q, hA, hB, hC, req::work_estimation, descr,
                        size_buf, buffer1, {}).wait();

    mkl::sparse::matmat(q, hA, hB, hC, req::get_compute_buf_size, descr,
                        size_buf, nullptr, {}).wait();
    auto* buffer2 = (void*)The_Arena()->alloc(std::size_t(*size_buf));
    mkl::sparse::matmat(q, hA, hB, hC, req::compute, descr,
                        size_buf, buffer2, {}).wait();

    mkl::sparse::matmat(q, hA, hB, hC, req::get_nnz, descr,
                        size_buf, nullptr, {}).wait();
    Long const c_nnz = *size_buf;
    AMREX_ALWAYS_ASSERT(c_nnz < Long(std::numeric_limits<I>::max()));

    C.mat.resize(c_nnz);
    C.col_index.resize(c_nnz);
    C.nnz = c_nnz;
#if defined(INTEL_MKL_VERSION) && (INTEL_MKL_VERSION < 20250300)
    mkl::sparse::set_csr_data(q, hC, I(nrows), I(ncols), mkl::index_base::zero,
                              C.row_offset.data(), C.col_index.data(), C.mat.data());
#else
    mkl::sparse::set_csr_data(q, hC, I(nrows), I(ncols), I(c_nnz), mkl::index_base::zero,
                              C.row_offset.data(), C.col_index.data(), C.mat.data());
#endif

    mkl::sparse::matmat(q, hA, hB, hC, req::finalize, descr,
                        size_buf, nullptr, {}).wait();

    mkl::sparse::release_matmat_descr(&descr);
    mkl::sparse::release_matrix_handle(q, &hA);
    mkl::sparse::release_matrix_handle(q, &hB);
    auto ev = mkl::sparse::release_matrix_handle(q, &hC);
    ev.wait();
    Gpu::streamSynchronize();
    The_Arena()->free(buffer1);
    The_Arena()->free(buffer2);
    The_Pinned_Arena()->free(size_buf);

    C.sort(); // oneMKL does not sort the output.

    return C;
}

#endif

/**
 * \brief Local CSR product C = A * B.
 *
 * `A` must have `nrows` rows and column indices in `[0, B.nrows)`; `B` must
 * have column indices in `[0, ncols)`. Rows of the result are sorted.
 */
template <typename T, template<typename> class V, typename I, typename F = NoRowPost>
CSR<T,V,I> spgemm_local (Long nrows, Long ncols,
                         CsrView<T const,I> const& A, CsrView<T const,I> const& B,
                         F const& row_post = {})
{
    AMREX_ASSERT(A.nrows == nrows);

    if (nrows <= 0 || ncols <= 0 || A.nnz <= 0 || B.nnz <= 0 || B.nrows <= 0) {
        return spgemm_empty<T,V,I>(nrows);
    }

#if !defined(AMREX_USE_GPU)
    return spgemm_local_cpu<T,V,I>(nrows, ncols, A, B, row_post);
#else
    amrex::ignore_unused(row_post);
#endif
#if defined(AMREX_USE_CUDA)
    return spgemm_local_cusparse<T,V,I>(nrows, ncols, A, B);
#elif defined(AMREX_USE_HIP)
    return spgemm_local_rocsparse<T,V,I>(nrows, ncols, A, B);
#elif defined(AMREX_USE_SYCL)
    return spgemm_local_onemkl<T,V,I>(nrows, ncols, A, B);
#endif
}

#ifdef AMREX_USE_MPI

// Row i of the result is row i of X0 followed by row i of X1, with column
// indices transformed by map0 and map1. Both views must have full row offsets.
template <typename T, template<typename> class V, typename I, typename F0, typename F1>
CSR<T,V,int> concat_csr_cols (Long nrows, CsrView<T const,I> const& X0, CsrView<T const,I> const& X1,
                              F0 const& map0, F1 const& map1)
{
    AMREX_ALWAYS_ASSERT(X0.nnz + X1.nnz < Long(std::numeric_limits<int>::max()));
    CSR<T,V,int> C;
    C.resize(nrows, X0.nnz + X1.nnz);
    int*  AMREX_RESTRICT crow = C.row_offset.data();
    int*  AMREX_RESTRICT ccol = C.col_index.data();
    T*    AMREX_RESTRICT cmat = C.mat.data();
    ParallelForOMP(nrows+1, [=] AMREX_GPU_DEVICE (Long i)
    {
        crow[i] = X0.row_offset[i] + X1.row_offset[i];
        if (i < nrows) {
            Long p = crow[i];
            for (Long q = X0.row_offset[i]; q < X0.row_offset[i+1]; ++q) {
                ccol[p] = map0(X0.col_index[q]);
                cmat[p] = X0.mat[q];
                ++p;
            }
            for (Long q = X1.row_offset[i]; q < X1.row_offset[i+1]; ++q) {
                ccol[p] = map1(X1.col_index[q]);
                cmat[p] = X1.mat[q];
                ++p;
            }
        }
    });
    return C;
}

// Append rows with global column indices, sorted within each row, to Bh
// whose first rows are already in the compact column space. Each row is
// rotated so that the local block [c0,c1) comes first, which keeps it
// sorted in the compact space [local | remote_union].
template <typename T, template<typename> class V>
void append_ext_rows (CSR<T,V,int>& Bh, Long const* ext_row_offset, Long const* ext_col,
                      T const* ext_mat, Long n_ext, Long nnz_ext,
                      Long c0, Long c1, Long const* ru, Long nru)
{
    Long const nb = Bh.nrows();
    Long const nnz_b = Bh.nnz;
    AMREX_ALWAYS_ASSERT(nnz_b + nnz_ext < Long(std::numeric_limits<int>::max()));
    Bh.mat.resize(nnz_b + nnz_ext);
    Bh.col_index.resize(nnz_b + nnz_ext);
    Bh.row_offset.resize(nb + n_ext + 1);
    Bh.nnz = nnz_b + nnz_ext;
    int*  AMREX_RESTRICT brow = Bh.row_offset.data();
    int*  AMREX_RESTRICT bcol = Bh.col_index.data();
    T*    AMREX_RESTRICT bmat = Bh.mat.data();
    Long const nlocal = c1 - c0;
    ParallelForOMP(n_ext, [=] AMREX_GPU_DEVICE (Long r)
    {
        Long const b = ext_row_offset[r];
        Long const e = ext_row_offset[r+1];
        brow[nb+r+1] = int(nnz_b + e);
        Long const p0 = amrex::lower_bound(ext_col+b, ext_col+e, c0) - ext_col;
        Long const p1 = amrex::lower_bound(ext_col+b, ext_col+e, c1) - ext_col;
        Long p = nnz_b + b;
        for (Long q = p0; q < p1; ++q) {
            bcol[p] = int(ext_col[q] - c0);
            bmat[p] = ext_mat[q];
            ++p;
        }
        for (Long q = b; q < e; ++q) {
            if (q < p0 || q >= p1) {
                bcol[p] = int(nlocal + (amrex::lower_bound(ru, ru+nru, ext_col[q]) - ru));
                bmat[p] = ext_mat[q];
                ++p;
            }
        }
    });
}

#endif

}

namespace amrex {

/**
 * \brief Sparse matrix-matrix multiplication, C = A * B.
 *
 * The row partition of `C` is that of `A`. The column partition of `A` is
 * set to the row partition of `B`, and the column partition of `B` to
 * `col_partition`, so each matrix can only ever be used with one column
 * partition. If `A` and `B` are the same object, `col_partition` must equal
 * its row partition.
 *
 * \param A First matrix.
 * \param B Second matrix.
 * \param col_partition Column partition of the product.
 *
 * \return Product of matrices `A` and `B`
 */
template <typename T, template <typename> class Allocator, typename F>
SpMatrix<T,Allocator>
SpGEMM (SpMatrix<T,Allocator> const& A, SpMatrix<T,Allocator> const& B,
        AlgPartition const& col_partition, F const& row_post);

template <typename T, template <typename> class Allocator>
SpMatrix<T,Allocator>
SpGEMM (SpMatrix<T,Allocator> const& A, SpMatrix<T,Allocator> const& B,
        AlgPartition const& col_partition)
{
    return SpGEMM(A, B, col_partition, detail::NoRowPost{});
}

/**
 * \brief SpGEMM with a callback on each row of the product.
 *
 * `row_post(i, col, val, n)` gets local row i with its n entries sorted by
 * column (compact column indices that sort like the global ones), may
 * change them in place, and returns the new number of entries, which must
 * stay sorted. CPU builds only; GPU builds accept only detail::NoRowPost.
 */
template <typename T, template <typename> class Allocator, typename F>
SpMatrix<T,Allocator>
SpGEMM (SpMatrix<T,Allocator> const& A, SpMatrix<T,Allocator> const& B,
        AlgPartition const& col_partition, F const& row_post)
{
    BL_PROFILE("SpGEMM");
#ifdef AMREX_USE_GPU
    static_assert(std::is_same_v<F,detail::NoRowPost>, "SpGEMM: row_post is CPU only");
    amrex::ignore_unused(row_post);
#endif

    using SpMat = SpMatrix<T,Allocator>;

    auto& Am = const_cast<SpMat&>(A);
    auto& Bm = const_cast<SpMat&>(B);
    Long const nrows = Am.numLocalRows();

#ifdef AMREX_USE_MPI
    using LongVec = typename SpMat::template container_type<Long>;

    Am.setColumnPartition(Bm.partition());
    Bm.setColumnPartition(col_partition);

    Long const nb = Bm.numLocalRows();
    Long const c0 = col_partition.globalRowBegin();
    Long const c1 = col_partition.globalRowEnd();
    Long const nlocal = c1 - c0;

    // Rows of B needed by A's off-diagonal columns, global column indices.
    typename SpMat::RemoteRowsMM ext;
    if (! detail::spmat_comm_is_local(Am.partition(), Bm.partition())) {
        ext = Am.fetch_remote_rows_mm(Bm);
    }

    if (ext.nrows == 0 && Am.m_remote_cols_v.empty() && Bm.m_remote_cols_v.empty()) {
        auto Ch = detail::spgemm_local<T,SpMat::template container_type,int>
            (nrows, nlocal, Am.m_csr_local.const_view(), Bm.m_csr_local.const_view(), row_post);
        SpMat C;
        C.define_split(Am.partition(), col_partition, std::move(Ch), nlocal, nullptr, 0);
        return C;
    }

    // Compact column space: [local columns | sorted remote columns].
    // TODO: this serial host loop over ext.nnz limits scaling on GPUs.
    Vector<Long> ru_h = Bm.m_remote_cols_v;
    for (Long i = 0; i < ext.nnz; ++i) {
        auto g = ext.col_index[i];
        if (g < c0 || g >= c1) { ru_h.push_back(g); }
    }
    RemoveDuplicates(ru_h);
    Long const nru = Long(ru_h.size());
    LongVec ru_d(nru);
    Gpu::copyAsync(Gpu::hostToDevice, ru_h.begin(), ru_h.end(), ru_d.begin());
    Long const* ru = ru_d.data();
    Long const ncols_hat = nlocal + nru;
    AMREX_ALWAYS_ASSERT(ncols_hat < Long(std::numeric_limits<int>::max()) &&
                        nb + Long(Am.m_remote_cols_v.size()) < Long(std::numeric_limits<int>::max()));

    using local_csr_type = typename SpMat::local_csr_type;
#ifndef AMREX_USE_GPU
    {
        // The rows of B fetched from other processes, in the compact
        // column space; B's own remote columns are mapped once each.
        local_csr_type E;
        E.resize(0, 0);
        E.row_offset[0] = 0;
        if (ext.nrows > 0) {
            detail::append_ext_rows(E, ext.row_offset.data(), ext.col_index, ext.mat,
                                    ext.nrows, ext.nnz, c0, c1, ru, nru);
        }
        ext.clear();
        Vector<int> b1map(Bm.m_remote_cols_v.size());
        for (Long j = 0; j < Long(b1map.size()); ++j) {
            b1map[j] = int(nlocal + (std::lower_bound(ru, ru+nru, Bm.m_remote_cols_v[j]) - ru));
        }
        auto const a1 = Am.m_remote_cols_v.empty() ? CsrView<T const,int>{}
                                                   : Am.remote_full_const_view();
        auto const b1 = Bm.m_remote_cols_v.empty() ? CsrView<T const,int>{}
                                                   : Bm.remote_full_const_view();
        auto Ch = detail::spgemm_split_cpu<T,SpMat::template container_type,int>
            (nrows, ncols_hat, nb, Am.m_csr_local.const_view(), a1,
             Bm.m_csr_local.const_view(), b1, b1map.data(), E.const_view(), row_post);
        SpMat C;
        C.define_split(Am.partition(), col_partition, std::move(Ch), nlocal, ru, nru);
        return C;
    }
#endif
    local_csr_type Ah = detail::concat_csr_cols<T,SpMat::template container_type>
        (nrows, Am.m_csr_local.const_view(), Am.remote_full_const_view(),
         [=] AMREX_GPU_DEVICE (int c) { return c; },
         [=] AMREX_GPU_DEVICE (int j) { return int(nb + j); });

#ifdef AMREX_USE_GPU
    Long const* b_rcols = Bm.m_remote_cols_dv.data();
#else
    Long const* b_rcols = Bm.m_remote_cols_v.data();
#endif
    local_csr_type Bh = detail::concat_csr_cols<T,SpMat::template container_type>
        (nb, Bm.m_csr_local.const_view(), Bm.remote_full_const_view(),
         [=] AMREX_GPU_DEVICE (int c) { return c; },
         [=] AMREX_GPU_DEVICE (int j) {
             return int(nlocal + (amrex::lower_bound(ru, ru+nru, b_rcols[j]) - ru)); });

    if (ext.nrows > 0) {
        LongVec ext_col_d(ext.nnz);
        Gpu::copyAsync(Gpu::hostToDevice, ext.col_index, ext.col_index+ext.nnz,
                       ext_col_d.begin());
        detail::append_ext_rows(Bh, ext.row_offset.data(), ext_col_d.data(), ext.mat,
                                ext.nrows, ext.nnz, c0, c1, ru, nru);
        Gpu::streamSynchronize();
    }
    ext.clear();

    local_csr_type Ch = detail::spgemm_local<T,SpMat::template container_type,int>
        (nrows, ncols_hat, Ah.const_view(), Bh.const_view());
    Ah = local_csr_type{};
    Bh = local_csr_type{};

    // The compact column space is [local | sorted remote], so the result
    // can be defined directly in split form.
    SpMat C;
    C.define_split(Am.partition(), col_partition, std::move(Ch), nlocal, ru, nru);
    return C;

#else

    Am.setColumnPartition(Bm.partition());
    Bm.setColumnPartition(col_partition);

    // One rank: local column indices are the global ones.
    Long const ncols = col_partition.numGlobalRows();
    AMREX_ALWAYS_ASSERT(ncols < Long(std::numeric_limits<int>::max()));
    auto Ch = detail::spgemm_local<T,SpMat::template container_type,int>
        (nrows, ncols, Am.m_csr_local.const_view(), Bm.m_csr_local.const_view(), row_post);
    SpMat C;
    C.define_split(Am.partition(), col_partition, std::move(Ch), ncols, nullptr, 0);
    return C;

#endif
}

/**
 * \brief Galerkin product R (A P), with the result of
 * SpGEMM(R, SpGEMM(A, P, col_partition), col_partition).
 *
 * In CPU builds A P is not stored in full: only its rows that other
 * processes need are formed as a matrix, and the local rows are computed
 * in chunks as the rows of R need them. The row partition of `R` must be
 * `col_partition`, and its column partition is set to the row partition of
 * `A`, as in SpGEMM.
 */
template <typename T, template <typename> class Allocator>
SpMatrix<T,Allocator>
RAP (SpMatrix<T,Allocator> const& R, SpMatrix<T,Allocator> const& A,
     SpMatrix<T,Allocator> const& P, AlgPartition const& col_partition)
{
    BL_PROFILE("RAP");
#ifdef AMREX_USE_GPU
    auto AP = SpGEMM(A, P, col_partition);
    return SpGEMM(R, AP, col_partition);
#else
    using SpMat = SpMatrix<T,Allocator>;
    using local_csr_type = typename SpMat::local_csr_type;
    auto& Rm = const_cast<SpMat&>(R);
    auto& Am = const_cast<SpMat&>(A);
    auto& Pm = const_cast<SpMat&>(P);
    Am.setColumnPartition(Pm.partition());
    Pm.setColumnPartition(col_partition);
    Rm.setColumnPartition(Am.partition());

    Long const nf = Am.numLocalRows();
    Long const nc = col_partition.numLocalRows();
    AMREX_ALWAYS_ASSERT(Rm.numLocalRows() == nc);

    detail::APRowsCpu<T,int> ap;
    ap.A0 = Am.m_csr_local.const_view();
    ap.P0 = Pm.m_csr_local.const_view();
    CsrView<T const,int> r1{};
    local_csr_type pe, ape;
    pe.resize(0, 0);
    pe.row_offset[0] = 0;
    ape.resize(0, 0);
    ape.row_offset[0] = 0;
    ap.PE = pe.const_view();
    Vector<int> p1map;
    Vector<Long> ru;

#ifdef AMREX_USE_MPI
    Long const c0 = col_partition.globalRowBegin();
    Long const c1 = col_partition.globalRowEnd();
    if (! Am.m_remote_cols_v.empty()) { ap.A1 = Am.remote_full_const_view(); }
    if (! Pm.m_remote_cols_v.empty()) { ap.P1 = Pm.remote_full_const_view(); }

    // Remote columns [ru_first, ru_last) as sorted unique global indices,
    // P's remote columns mapped into the compact space of `cols`, and the
    // fetched rows `ext` in that space.
    auto add_remote = [&] (Vector<Long>& cols, typename SpMat::RemoteRowsMM const& ext) {
        for (Long i = 0; i < ext.nnz; ++i) {
            auto g = ext.col_index[i];
            if (g < c0 || g >= c1) { cols.push_back(g); }
        }
        RemoveDuplicates(cols);
    };
    auto compact_ext = [&] (local_csr_type& E, typename SpMat::RemoteRowsMM const& ext,
                            Vector<Long> const& cols) {
        E.resize(0, 0);
        E.row_offset[0] = 0;
        if (ext.nrows > 0) {
            detail::append_ext_rows(E, ext.row_offset.data(), ext.col_index, ext.mat,
                                    ext.nrows, ext.nnz, c0, c1, cols.data(), Long(cols.size()));
        }
    };
    auto map_p1 = [&] (Vector<Long> const& cols) {
        p1map.resize(Pm.m_remote_cols_v.size());
        for (Long j = 0; j < Long(p1map.size()); ++j) {
            p1map[j] = int(nc + (std::lower_bound(cols.begin(), cols.end(),
                                                  Pm.m_remote_cols_v[j]) - cols.begin()));
        }
        ap.p1map = p1map.data();
    };

    // Rows of P needed by A's remote columns.
    typename SpMat::RemoteRowsMM pext;
    if (! detail::spmat_comm_is_local(Am.partition(), Pm.partition())) {
        pext = Am.fetch_remote_rows_mm(Pm);
    }
    Vector<Long> ru_ap = Pm.m_remote_cols_v;
    add_remote(ru_ap, pext);

    // The rows of A P that other processes need: those where P has remote
    // columns, since the transpose of P then has a remote column there.
    SpMat APb;
    {
        compact_ext(pe, pext, ru_ap);
        ap.PE = pe.const_view();
        map_p1(ru_ap);
        Long const ncols = nc + Long(ru_ap.size());
        AMREX_ALWAYS_ASSERT(ncols < Long(std::numeric_limits<int>::max()));
        auto has_remote = [&] (Long i) {
            return ap.P1.nnz > 0 && ap.P1.row_offset[i+1] > ap.P1.row_offset[i];
        };
        auto csr = detail::csr_from_rows_cpu<T,SpMat::template container_type,int>
            (nf, [&] (Long i) { return has_remote(i) ? ap.maxlen(i) : Long(0); },
             [&] (int) {
            return [&, marker = Vector<int>(ncols, -1), col = Vector<int>(), val = Vector<T>(),
                    tmp = Vector<std::pair<int,T>>()]
                (Long i, int const*& rc, T const*& rv) mutable -> Long
            {
                Long n = 0;
                if (has_remote(i)) {
                    Long const maxlen = ap.maxlen(i);
                    if (Long(col.size()) < maxlen) {
                        col.resize(maxlen);
                        val.resize(maxlen);
                    }
                    n = ap.row(i, marker.data(), col.data(), val.data(), 0);
                    // fetch_remote_rows_mm sends rows sorted by column.
                    detail::sort_row_cpu(col.data(), val.data(), n, tmp);
                }
                rc = col.data();
                rv = val.data();
                return n;
            };
        });
        APb.define_split(Am.partition(), col_partition, std::move(csr), nc,
                         ru_ap.data(), Long(ru_ap.size()));
    }

    // Rows of A P needed by R's remote columns.
    typename SpMat::RemoteRowsMM apext;
    if (! detail::spmat_comm_is_local(Rm.partition(), APb.partition())) {
        apext = Rm.fetch_remote_rows_mm(APb);
    }
    ru = ru_ap;
    add_remote(ru, apext);
    AMREX_ALWAYS_ASSERT(nc + Long(ru.size()) < Long(std::numeric_limits<int>::max()));
    compact_ext(pe, pext, ru);
    compact_ext(ape, apext, ru);
    ap.PE = pe.const_view();
    map_p1(ru);
    pext.clear();
    apext.clear();
    if (! Rm.m_remote_cols_v.empty()) { r1 = Rm.remote_full_const_view(); }
#endif

    Long const nru = Long(ru.size());
    auto csr = detail::rap_local_cpu<T,SpMat::template container_type,int>
        (nf, nc, nc+nru, Rm.m_csr_local.const_view(), r1, ape.const_view(), ap);
    SpMat C;
    C.define_split(Rm.partition(), col_partition, std::move(csr), nc, ru.data(), nru);
    return C;
#endif
}

}

#endif
