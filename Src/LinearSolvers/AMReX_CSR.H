#ifndef AMREX_CSR_H_
#define AMREX_CSR_H_
#include <AMReX_Config.H>

#include <AMReX_Gpu.H>
#include <AMReX_INT.H>
#include <AMReX_OpenMP.H>

#if defined(AMREX_USE_CUDA)
#include <cub/cub.cuh> // for Clang
#endif

#include <algorithm>
#include <climits>
#include <type_traits>

namespace amrex {

template <typename T>
struct CsrView {
    using U = std::conditional_t<std::is_const_v<T>, Long const, Long>;
    T* AMREX_RESTRICT mat = nullptr;
    U* AMREX_RESTRICT col_index = nullptr;
    U* AMREX_RESTRICT row_offset = nullptr;
    Long nnz = 0;
    Long nrows = 0;
};

template <typename T, template <typename> class V>
struct CSR {
    V<T> mat;
    V<Long> col_index;
    V<Long> row_offset;
    Long nnz = 0;

    [[nodiscard]] Long nrows () const {
        return row_offset.empty() ? Long(0) : Long(row_offset.size())-1;
    }

    void resize (Long num_rows, Long num_non_zeros) {
        mat.resize(num_non_zeros);
        col_index.resize(num_non_zeros);
        row_offset.resize(num_rows+1);
        nnz = num_non_zeros;
    }

    CsrView<T> view () {
        return CsrView<T>{mat.data(), col_index.data(), row_offset.data(),
                          nnz, Long(row_offset.size())-1};
    }

    [[nodiscard]] CsrView<T const> view () const {
        return CsrView<T>{mat.data(), col_index.data(), row_offset.data(),
                          nnz, Long(row_offset.size())-1};
    }

    [[nodiscard]] CsrView<T const> const_view () const {
        return CsrView<T const>{mat.data(), col_index.data(), row_offset.data(),
                                nnz, Long(row_offset.size())-1};
    }

    void sort ();

    void sort_on_host ();
};

template <typename C, typename T, template<typename> class AD, template<typename> class AS,
          std::enable_if_t<std::is_same_v<C,Gpu::HostToDevice> ||
                           std::is_same_v<C,Gpu::DeviceToHost> ||
                           std::is_same_v<C,Gpu::DeviceToDevice>, int> = 0>
void duplicateCSR (C c, CSR<T,AD>& dst, CSR<T,AS> const& src)
{
    dst.mat.resize(src.mat.size());
    dst.col_index.resize(src.col_index.size());
    dst.row_offset.resize(src.row_offset.size());
    Gpu::copyAsync(c,
                   src.mat.begin(),
                   src.mat.end(),
                   dst.mat.begin());
    Gpu::copyAsync(c,
                   src.col_index.begin(),
                   src.col_index.end(),
                   dst.col_index.begin());
    Gpu::copyAsync(c,
                   src.row_offset.begin(),
                   src.row_offset.end(),
                   dst.row_offset.begin());
    dst.nnz = src.nnz;
}

template <typename T, template <typename> class V>
void CSR<T,V>::sort ()
{
    if (nnz <= 0) { return; }

#ifdef AMREX_USE_GPU

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

    // The function is synchronous. If that is no longer the case, we might
    // need to update SpMatrix::define.

    constexpr int nthreads = 256;
    constexpr int nwarps_per_block = nthreads / Gpu::Device::warp_size;

    AMREX_ALWAYS_ASSERT((nrows()+nwarps_per_block-1) < Long(std::numeric_limits<int>::max()));

    auto nr = int(nrows());
    int nblocks = (nr + nwarps_per_block-1) / nwarps_per_block;
    auto const& stream = Gpu::gpuStream();

    auto* pmat = mat.data();
    auto* pcol = col_index.data();
    auto* prow = row_offset.data();

    Gpu::Buffer<int> needs_fallback({0});
    auto* d_needs_fallback = needs_fallback.data();

    amrex::launch_global<nthreads><<<nblocks, nthreads, 0, stream>>>
        ([=] AMREX_GPU_DEVICE () noexcept
    {
        int wid = int(threadIdx.x)/Gpu::Device::warp_size;
        int r = int(blockIdx.x)*nwarps_per_block + wid;
        if (r >= nr) return;

        Long const b = prow[r];
        Long const e = prow[r+1];
        auto const len = int(e - b);

        if (len <= 1) return;

        int lane = int(threadIdx.x) - wid * Gpu::Device::warp_size;

        bool sorted = true;
        for (Long i = lane + 1; i < len; i += Gpu::Device::warp_size) {
            sorted = sorted && (pcol[b+i-1] <= pcol[b+i]);
        }
#if defined(AMREX_USE_CUDA)
        if (__all_sync(0xffffffff, sorted)) { return; }
#else
        if (__all(sorted)) { return; }
#endif

        constexpr int ITEMS_PER_THREAD = AMREX_HIP_OR_CUDA(2,4);
        constexpr int ITEMS_PER_WARP = Gpu::Device::warp_size * ITEMS_PER_THREAD;

        if (len <= ITEMS_PER_WARP)
        {
#if defined(AMREX_USE_CUDA)
            using WarpSort = cub::WarpMergeSort<Long, ITEMS_PER_THREAD, Gpu::Device::warp_size, T>;
            __shared__ typename WarpSort::TempStorage temp_storage[nwarps_per_block];
#elif defined(AMREX_USE_HIP)
            using WarpSort = rocprim::warp_sort<Long, Gpu::Device::warp_size, T>;
            __shared__ typename WarpSort::storage_type temp_storage[nwarps_per_block];
#endif

            Long keys[ITEMS_PER_THREAD];
            T values[ITEMS_PER_THREAD];

            #pragma unroll
            for (int i = 0; i < ITEMS_PER_THREAD; ++i) {
                int idx = lane * ITEMS_PER_THREAD + i;
                if (idx < len) {
                    keys[i] = pcol[b + idx];
                    values[i] = pmat[b + idx];
                } else {
                    keys[i] = std::numeric_limits<Long>::max();
                    values[i] = T{};
                }
            }

            AMREX_HIP_OR_CUDA(
                WarpSort{}.sort(keys, values, temp_storage[wid]),
                WarpSort(temp_storage[wid]).Sort(
                    keys, values, [](Long x, Long y) {return x < y;}));

            #pragma unroll
            for (int i = 0; i < ITEMS_PER_THREAD; ++i) {
                int idx = lane * ITEMS_PER_THREAD + i;
                if (idx < len) {
                    pcol[b + idx] = keys[i];
                    pmat[b + idx] = values[i];
                }
            }
        } else {
            if (lane == 0) {
                Gpu::Atomic::AddNoRet(d_needs_fallback, 1);
            }
        }
    });

    auto* h_needs_fallback = needs_fallback.copyToHost();

    if (*h_needs_fallback)
    {
        V<Long> col_index_out(col_index.size());
        V<T> mat_out(mat.size());
        auto* d_col_out = col_index_out.data();
        auto* d_val_out = mat_out.data();

        std::size_t temp_bytes = 0;

        AMREX_GPU_SAFE_CALL(AMREX_HIP_OR_CUDA(
                                rocprim::segmented_radix_sort_pairs,
                                cub::DeviceSegmentedRadixSort::SortPairs)
                            (nullptr, temp_bytes, pcol, d_col_out, pmat, d_val_out,
                             nnz, nr, prow, prow+1, 0, int(sizeof(Long)*CHAR_BIT),
                             stream));

        auto* d_temp = (void*) The_Arena()->alloc(temp_bytes);

        AMREX_GPU_SAFE_CALL(AMREX_HIP_OR_CUDA(
                                rocprim::segmented_radix_sort_pairs,
                                cub::DeviceSegmentedRadixSort::SortPairs)
                            (d_temp, temp_bytes, pcol, d_col_out, pmat, d_val_out,
                             nnz, nr, prow, prow+1, 0, int(sizeof(Long)*CHAR_BIT),
                             stream));

        std::swap(col_index, col_index_out);
        std::swap(mat, mat_out);

        Gpu::streamSynchronize();
        The_Arena()->free(d_temp);
    }

    // let's test both by print matrix out to see if it's sorted.

    AMREX_GPU_ERROR_CHECK();

#elif defined(AMREX_USE_SYCL)

    // xxxxx TODO SYCL: Let's not worry about performance for now.
    CSR<T,Gpu::PinnedVector> h_csr;
    duplicateCSR(Gpu::deviceToHost, h_csr, *this);
    Gpu::streamSynchronize();
    h_csr.sort_on_host();
    duplicateCSR(Gpu::hostToDevice, *this, h_csr);
    Gpu::streamSynchronize();

#endif

#else

    sort_on_host();

#endif
}

template <typename T, template <typename> class V>
void CSR<T,V>::sort_on_host ()
{
    if (nnz <= 0) { return; }

    constexpr int SMALL = 128;

    Long nr = nrows();

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    {
        V<Long> lcols;
        V<T   > lvals;
        V<int > perm;

        Long scols[SMALL];
        T    svals[SMALL];

#ifdef AMREX_USE_OMP
#pragma omp for
#endif
        for (Long r = 0; r < nr; ++r) {
            Long const b = row_offset[r  ];
            Long const e = row_offset[r+1];
            auto const len = int(e - b);

            if (len <= 1) { continue; }

            bool sorted = true;
            for (int i = 1; i < len; ++i) {
                if (col_index[b+i-1] > col_index[b+i]) {
                    sorted = false;
                    break;
                }
            }
            if (sorted) { continue; }

            if (len <= SMALL) {
                // Insertion sort using arrays on stack
                for (int i = 0; i < len; ++i) {
                    scols[i] = col_index[b+i];
                    svals[i] = mat      [b+i];
                }
                for (int i = 1; i < len; ++i) {
                    auto c = scols[i];
                    auto v = svals[i];
                    auto j = i;
                    while (j > 0 && scols[j-1] > c) {
                        scols[j] = scols[j-1];
                        svals[j] = svals[j-1];
                        --j;
                    }
                    scols[j] = c;
                    svals[j] = v;
                }
                for (int i = 0; i < len; ++i) {
                    col_index[b+i] = scols[i];
                    mat      [b+i] = svals[i];
                }
            } else {
                lcols.resize(len);
                lvals.resize(len);
                perm.resize(len);

                for (int i = 0; i < len; ++i) {
                    lcols[i] = col_index[b+i];
                    lvals[i] = mat      [b+i];
                    perm [i] = i;
                }

                std::sort(perm.begin(), perm.end(),
                          [&] (int i0, int i1) {
                              return lcols[i0] < lcols[i1];
                          });

                for (int out = 0; out < len; ++out) {
                    auto const in = perm[out];
                    col_index[b+out] = lcols[in];
                    mat      [b+out] = lvals[in];
                }
            }
        }
    }
}

}

#endif
