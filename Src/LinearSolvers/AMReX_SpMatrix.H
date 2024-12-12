#ifndef AMREX_SP_MATRIX_H_
#define AMREX_SP_MATRIX_H_
#include <AMReX_Config.H>

#include <AMReX_AlgPartition.H>
#include <AMReX_AlgVector.H>
#include <AMReX_Gpu.H>
#include <AMReX_INT.H>
#include <AMReX_Scan.H>

#include <fstream>
#include <string>
#include <type_traits>

namespace amrex {

template <typename T, template<typename> class Allocator = DefaultAllocator>
class SpMatrix
{
public:
    using value_type = T;

    SpMatrix () = default;

    SpMatrix (AlgPartition partition, int nnz);

    SpMatrix (SpMatrix const&) = delete;
    SpMatrix& operator= (SpMatrix const&) = delete;

    SpMatrix (SpMatrix &&) = default;
    SpMatrix& operator= (SpMatrix &&) = default;

    ~SpMatrix () = default;

    void define (AlgPartition partition, int nnz);

    [[nodiscard]] AlgPartition const& partition () const { return m_partition; }

    [[nodiscard]] Long numLocalRows () const { return m_row_end - m_row_begin; }
    [[nodiscard]] Long numGlobalRows () const { return m_partition.numGlobalRows(); }
    [[nodiscard]] Long numLocalNonZero () const { return m_data.nnz; }

    //! Inclusive global index begin.
    [[nodiscard]] Long globalRowBegin () const { return m_row_begin; }
    //! Exclusive global index end.
    [[nodiscard]] Long globalRowEnd () const { return m_row_end; }

    [[nodiscard]] T const* data () const { return m_data.mat.data(); }
    [[nodiscard]] T      * data ()       { return m_data.mat.data(); }
    [[nodiscard]] Long const* columnIndex () const { return m_data.col_index.data(); }
    [[nodiscard]] Long      * columnIndex ()       { return m_data.col_index.data(); }
    [[nodiscard]] Long const* rowOffset () const { return m_data.row_offset.data(); }
    [[nodiscard]] Long      * rowOffset ()       { return m_data.row_offset.data(); }

    void printToFile (std::string const& file) const;

    template <typename F>
    void setVal (F const& f);

    [[nodiscard]] AlgVector<T> const& diagonalVector () const;

    template <typename U> friend void SpMV(AlgVector<U>& y, SpMatrix<U> const& A, AlgVector<U> const& x);

    //! Private function, but public for cuda
    void define_doit (int nnz);

#ifdef AMREX_USE_MPI
    //! Private function, but public for cuda
    void prepare_comm ();
    void pack_buffer (AlgVector<T> const& v);
    void unpack_buffer (AlgVector<T>& v);
#endif

private:

    void startComm (AlgVector<T> const& x);
    void finishComm (AlgVector<T>& y);

    template <class U> using DVec = PODVector<U,Allocator<U> >;

    template <template <typename> class V>
    struct CSR {
        V<T> mat;
        V<Long> col_index;
        V<Long> row_offset;
        Long nnz = -1;
    };

    AlgPartition m_partition;
    Long m_row_begin = 0;
    Long m_row_end = 0;
    CSR<DVec> m_data;

    mutable AlgVector<T> m_diagonal;

#ifdef AMREX_USE_MPI
    CSR<DVec> m_data_remote;

#ifdef AMREX_USE_GPU
    Gpu::PinnedVector<Long> m_remote_cols;
#endif

    DVec<Long> m_rtol;

    Vector<int> m_send_to;
    Vector<int> m_send_counts;
    Gpu::DeviceVector<Long> m_send_indices;

    Vector<int> m_recv_from;
    Vector<int> m_recv_counts;

    Vector<MPI_Request> m_send_reqs;
    T* m_send_buffer = nullptr;
    Long m_total_counts_send = 0;

    Vector<MPI_Request> m_recv_reqs;
    T* m_recv_buffer = nullptr;
    Long m_total_counts_recv = 0;

    bool m_comm_prepared = false;
#endif

    bool m_shifted = false;
};

template <typename T, template<typename> class Allocator>
SpMatrix<T,Allocator>::SpMatrix (AlgPartition partition, int nnz)
    : m_partition(std::move(partition)),
      m_row_begin(m_partition[ParallelDescriptor::MyProc()]),
      m_row_end(m_partition[ParallelDescriptor::MyProc()+1])
{
    static_assert(std::is_floating_point<T>::value, "SpMatrix is for floating point type only");
    define_doit(nnz);
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::define (AlgPartition partition, int nnz)
{
    m_partition = std::move(partition);
    m_row_begin = m_partition[ParallelDescriptor::MyProc()];
    m_row_end = m_partition[ParallelDescriptor::MyProc()+1];
    define_doit(nnz);
}

template <typename T, template<typename> class Allocator>
void
SpMatrix<T,Allocator>::define_doit (int nnz)
{
    Long nlocalrows = this->numLocalRows();
    Long total_nnz = nlocalrows*nnz;
    m_data.mat.resize(total_nnz);
    m_data.col_index.resize(total_nnz);
    m_data.row_offset.resize(nlocalrows+1);
    m_data.nnz = total_nnz;

    auto* poffset = m_data.row_offset.data();
    ParallelFor(nlocalrows+1, [=] AMREX_GPU_DEVICE (Long lrow) noexcept
    {
        poffset[lrow] = lrow*nnz;
    });
}

template <typename T, template<typename> class Allocator>
void
SpMatrix<T,Allocator>::printToFile (std::string const& file) const
{
    // xxxxx TODO: This function only prints square part of the local rows,
    // not the full rows.

#ifdef AMREX_USE_GPU
    CSR<Gpu::PinnedVector> csr;
    csr.mat.resize(m_data.mat.size());
    csr.col_index.resize(m_data.col_index.size());
    csr.row_offset.resize(m_data.row_offset.size());
    Gpu::copyAsync(Gpu::deviceToHost, m_data.mat.begin(), m_data.mat.end(), csr.mat.begin());
    Gpu::copyAsync(Gpu::deviceToHost, m_data.col_index.begin(), m_data.col_index.end(), csr.col_index.begin());
    Gpu::copyAsync(Gpu::deviceToHost, m_data.row_offset.begin(), m_data.row_offset.end(), csr.row_offset.begin());
    csr.nnz = m_data.nnz;
    Gpu::streamSynchronize();
#else
    auto const& csr = m_data;
#endif

    const Long lrow_begin = m_partition[ParallelDescriptor::MyProc()];
    std::ofstream ofs(file+"."+std::to_string(ParallelDescriptor::MyProc()));
    ofs << m_row_begin << " " << m_row_end << " " << csr.nnz << "\n";
    for (Long i = 0, nrows = numLocalRows(); i < nrows; ++i) {
        Long nnz_row = csr.row_offset[i+1] - csr.row_offset[i];
        T    const* mat = csr.mat.data()       + csr.row_offset[i];
        Long const* col = csr.col_index.data() + csr.row_offset[i];
        for (Long j = 0; j < nnz_row; ++j) {
            ofs << i+lrow_begin << " " << col[j] << " " << mat[j] << "\n";
        }
    }
}

template <typename T, template<typename> class Allocator>
template <typename F>
void SpMatrix<T,Allocator>::setVal (F const& f)
{
    // xxxxx TODO: We can try to optimize this later by using shared memory.

    Long nlocalrows = this->numLocalRows();
    Long rowbegin = this->globalRowBegin();
    auto* pmat = m_data.mat.data();
    auto* pcolindex = m_data.col_index.data();
    auto* prowoffset = m_data.row_offset.data();
    ParallelFor(nlocalrows, [=] AMREX_GPU_DEVICE (int lrow) noexcept
    {
        f(rowbegin+lrow, pcolindex+prowoffset[lrow], pmat+prowoffset[lrow]);
    });
}

template <typename T, template<typename> class Allocator>
AlgVector<T> const& SpMatrix<T,Allocator>::diagonalVector () const
{
    if (m_diagonal.empty()) {
        m_diagonal.define(this->partition());
        auto* AMREX_RESTRICT p = m_diagonal.data();
        auto const* AMREX_RESTRICT mat = m_data.mat.data();
        auto const* AMREX_RESTRICT col = m_data.col_index.data();
        auto const* AMREX_RESTRICT row = m_data.row_offset.data();
        auto offset = m_shifted ? Long(0) : m_row_begin;
        Long nrows = this->numLocalRows();
        ParallelFor(nrows, [=] AMREX_GPU_DEVICE (Long i)
        {
            T d = 0;
            for (Long j = row[i]; j < row[i+1]; ++j) {
                if (i == col[j] - offset) {
                    d = mat[j];
                    break;
                }
            }
            p[i] = d;
        });
    }
    return m_diagonal;
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::startComm (AlgVector<T> const& x)
{
#ifndef AMREX_USE_MPI
    amrex::ignore_unused(x);
#else
    if (this->partition().numActiveProcs() <= 1) { return; }

    this->prepare_comm();

    auto const mpi_tag = ParallelDescriptor::SeqNum();
    auto const mpi_t_type = ParallelDescriptor::Mpi_typemap<T>::type(); // NOLINT(readability-qualified-auto)
    auto const mpi_comm = ParallelContext::CommunicatorSub(); // NOLINT(readability-qualified-auto)

    auto const nrecvs = int(m_recv_from.size());
    if (nrecvs > 0) {
        m_recv_buffer = (T*)The_Comms_Arena()->alloc(sizeof(T)*m_total_counts_recv);
        m_recv_reqs.resize(nrecvs, MPI_REQUEST_NULL);
        auto* p_recv = m_recv_buffer;
        for (int irecv = 0; irecv < nrecvs; ++irecv) {
            BL_MPI_REQUIRE(MPI_Irecv(p_recv,
                                     m_recv_counts[irecv], mpi_t_type,
                                     m_recv_from[irecv], mpi_tag, mpi_comm,
                                     &(m_recv_reqs[irecv])));
            p_recv += m_recv_counts[irecv];
        }
        AMREX_ASSERT(p_recv == m_recv_buffer + m_total_counts_recv);
    }

    auto const nsends = int(m_send_to.size());
    if (nsends > 0) {
        m_send_buffer = (T*)The_Comms_Arena()->alloc(sizeof(T)*m_total_counts_send);

        pack_buffer(x);
        Gpu::streamSynchronize();

        m_send_reqs.resize(nsends, MPI_REQUEST_NULL);
        auto* p_send = m_send_buffer;
        for (int isend = 0; isend < nsends; ++isend) {
            auto count = m_send_counts[isend];
            BL_MPI_REQUIRE(MPI_Isend(p_send, count, mpi_t_type, m_send_to[isend],
                                     mpi_tag, mpi_comm, &(m_send_reqs[isend])));
            p_send += count;
        }
        AMREX_ASSERT(p_send == m_send_buffer + m_total_counts_send);
    }
#endif
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::finishComm (AlgVector<T>& y)
{
    if (this->numLocalRows() == 0) { return; }

#ifndef AMREX_USE_MPI
    amrex::ignore_unused(y);
#else
    if (this->partition().numActiveProcs() <= 1) { return; }

    if ( ! m_recv_reqs.empty()) {
        Vector<MPI_Status> mpi_statuses(m_recv_reqs.size());
        BL_MPI_REQUIRE(MPI_Waitall(int(m_recv_reqs.size()),
                                   m_recv_reqs.data(),
                                   mpi_statuses.data()));
    }

    unpack_buffer(y);

    if ( ! m_send_reqs.empty()) {
        Vector<MPI_Status> mpi_statuses(m_send_reqs.size());
        BL_MPI_REQUIRE(MPI_Waitall(int(m_send_reqs.size()),
                                   m_send_reqs.data(),
                                   mpi_statuses.data()));
    }

    Gpu::streamSynchronize();
    The_Comms_Arena()->free(m_send_buffer);
    The_Comms_Arena()->free(m_recv_buffer);
    m_send_reqs.clear();
    m_recv_reqs.clear();
#endif
}

#ifdef AMREX_USE_MPI

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::prepare_comm ()
{
    if (m_comm_prepared) { return; }

    // This function needs to be safe when nnz is zero.

    // xxxxx TODO: check there is no int overflow.

    const int nprocs = ParallelContext::NProcsSub();

    // First, we need to split the matrix into two parts, a square matrix
    // for pure local operations and another part for remote operations.

    Long all_nnz = m_data.nnz;
    Long local_nnz;
    Gpu::DeviceVector<Long> pfsum(all_nnz);
    auto* p_pfsum = pfsum.data();
    auto row_begin = m_row_begin;
    auto row_end = m_row_end;
    if (m_data.nnz < Long(std::numeric_limits<int>::max())) {
        auto const* pcol = m_data.col_index.data();
        local_nnz = Scan::PrefixSum<int>(int(all_nnz),
                                         [=] AMREX_GPU_DEVICE (int i) -> int {
                                             return (pcol[i] >= row_begin &&
                                                     pcol[i] <  row_end); },
                                         [=] AMREX_GPU_DEVICE (int i, int const& x) {
                                             p_pfsum[i] = x; },
                                         Scan::Type::exclusive, Scan::retSum);
    } else {
        auto const* pcol = m_data.col_index.data();
        local_nnz = Scan::PrefixSum<Long>(all_nnz,
                                          [=] AMREX_GPU_DEVICE (Long i) -> Long {
                                              return (pcol[i] >= row_begin &&
                                                      pcol[i] <  row_end); },
                                          [=] AMREX_GPU_DEVICE (Long i, Long const& x) {
                                              p_pfsum[i] = x; },
                                          Scan::Type::exclusive, Scan::retSum);
    }

    m_data.nnz = local_nnz;
    Long remote_nnz = all_nnz - local_nnz;
    m_data_remote.nnz = remote_nnz;

    Vector<Vector<Long>>unique_remote_cols_vv(nprocs);
    Vector<Long> unique_remote_cols_v;

    if (local_nnz != all_nnz) {
        m_data_remote.mat.resize(remote_nnz);
        m_data_remote.col_index.resize(remote_nnz);
        DVec<T> new_mat(local_nnz);
        DVec<Long> new_col(local_nnz);
        auto const* pmat = m_data.mat.data();
        auto const* pcol = m_data.col_index.data();
        auto* pmat_l = new_mat.data();
        auto* pcol_l = new_col.data();
        auto* pmat_r = m_data_remote.mat.data();
        auto* pcol_r = m_data_remote.col_index.data();
        ParallelFor(all_nnz, [=] AMREX_GPU_DEVICE (Long i)
        {
            auto ps = p_pfsum[i];
            auto local = (pcol[i] >= row_begin &&
                          pcol[i] <  row_end);
            if (local) {
                pmat_l[ps] = pmat[i];
                pcol_l[ps] = pcol[i] - row_begin; // shift the column index to local
            } else {
                pmat_r[i-ps] = pmat[i];
                pcol_r[i-ps] = pcol[i];
            }
        });
        m_shifted = true;
        auto noffset = Long(m_data.row_offset.size());
        auto* pro = m_data.row_offset.data();
        m_data_remote.row_offset.resize(noffset);
        auto* pro_r = m_data_remote.row_offset.data();
        ParallelFor(noffset, [=] AMREX_GPU_DEVICE (Long i)
        {
            if (i < noffset-1) {
                auto ro_l = p_pfsum[pro[i]];
                pro_r[i] = pro[i] - ro_l;
                pro[i] = ro_l;
            } else {
                pro[i] = local_nnz;
                pro_r[i] = remote_nnz;
            }
        });
        Gpu::streamSynchronize();
        m_data.mat.swap(new_mat);
        m_data.col_index.swap(new_col);

        // In the remote part, it's expected that some rows don't have
        // non-zeros. So we trim them off.
        {
            Long old_size = m_data_remote.row_offset.size();
            m_rtol.resize(old_size-1);
            auto* p_rtol = m_rtol.data();
            DVec<Long> trimmed_row_offset(old_size);
            auto const* p_ro = m_data_remote.row_offset.data();
            auto* p_tro = trimmed_row_offset.data();
            Long new_size;
            if (old_size < Long(std::numeric_limits<int>::max())) {
                // This is basically std::unique.
                new_size = Scan::PrefixSum<int>(int(old_size),
                                                [=] AMREX_GPU_DEVICE (int i) -> int {
                                                    if (i+1 < old_size) {
                                                        return (p_ro[i+1] > p_ro[i]);
                                                    } else {
                                                        return 1;
                                                    }
                                                },
                                                [=] AMREX_GPU_DEVICE (int i, int const& x) {
                                                    if (i == 0) {
                                                        p_tro[0] = 0;
                                                    } else if (p_ro[i] > p_ro[i-1]) {
                                                        p_tro[x] = p_ro[i];
                                                    }
                                                    if ((i+1 < old_size) &&
                                                        p_ro[i+1] > p_ro[i])
                                                    {
                                                        p_rtol[x] = i;
                                                    }
                                                },
                                                Scan::Type::exclusive, Scan::retSum);
            } else {
                // This is basically std::unique.
                new_size = Scan::PrefixSum<Long>(old_size,
                                                [=] AMREX_GPU_DEVICE (Long i) -> Long {
                                                    if (i+1 < old_size) {
                                                        return (p_ro[i+1] > p_ro[i]);
                                                    } else {
                                                        return 1;
                                                    }
                                                },
                                                [=] AMREX_GPU_DEVICE (Long i, Long const& x) {
                                                    if (i == 0) {
                                                        p_tro[0] = 0;
                                                    } else if (p_ro[i] > p_ro[i-1]) {
                                                        p_tro[x] = p_ro[i];
                                                    }
                                                    if ((i+1 < old_size) &&
                                                        p_ro[i+1] > p_ro[i])
                                                    {
                                                        p_rtol[x] = i;
                                                    }
                                                },
                                                Scan::Type::exclusive, Scan::retSum);
            }

            m_rtol.resize(new_size-1);
            trimmed_row_offset.resize(new_size);
#ifdef AMREX_USE_GPU
            m_rtol.shrink_to_fit();
            trimmed_row_offset.shrink_to_fit();
#endif
            m_data_remote.row_offset.swap(trimmed_row_offset);
        }

#ifdef AMREX_USE_GPU
        m_remote_cols.resize(m_data_remote.col_index.size());
        Gpu::copyAsync(Gpu::deviceToHost, m_data_remote.col_index.begin(),
                                          m_data_remote.col_index.end(),
                                          m_remote_cols.begin());
        Gpu::streamSynchronize();
#else
        auto const& m_remote_cols = m_data_remote.col_index;
#endif

        unique_remote_cols_v.resize(m_remote_cols.size());
        std::partial_sort_copy(m_remote_cols.begin(),
                               m_remote_cols.end(),
                               unique_remote_cols_v.begin(),
                               unique_remote_cols_v.end());
        amrex::RemoveDuplicates(unique_remote_cols_v);

        m_total_counts_recv = Long(unique_remote_cols_v.size());

        // Note that amrex::RemoveDuplicates sorts the data.
        auto const& rows = this->m_partition.dataVector();
        auto it = rows.cbegin();
        for (auto c : unique_remote_cols_v) {
            it = std::find_if(it, rows.cend(), [&] (auto x) { return x > c; });
            if (it != rows.cend()) {
                int iproc = int(std::distance(rows.cbegin(),it)) - 1;
                unique_remote_cols_vv[iproc].push_back(c);
            } else {
                amrex::Abort("SpMatrix::prepare_comm: how did this happen?");
            }
        }
    }

    // Need to make plans for MPI
    auto const mpi_tag = ParallelDescriptor::SeqNum();
    auto const mpi_int = ParallelDescriptor::Mpi_typemap<int>::type(); // NOLINT(readability-qualified-auto)
    auto const mpi_long = ParallelDescriptor::Mpi_typemap<Long>::type(); // NOLINT(readability-qualified-auto)
    auto const mpi_comm = ParallelContext::CommunicatorSub(); // NOLINT(readability-qualified-auto)

    amrex::Vector<int> need_from(nprocs);
    for (int iproc = 0; iproc < nprocs; ++iproc) {
        need_from[iproc] = unique_remote_cols_vv[iproc].empty() ? 0 : 1;
    }
    amrex::Vector<int> reduce_scatter_counts(nprocs,1);
    int nsends = 0;
    BL_MPI_REQUIRE(MPI_Reduce_scatter
                   (need_from.data(), &nsends, reduce_scatter_counts.data(),
                    mpi_int, MPI_SUM, mpi_comm));

    // nsends is the number of processes that need data from me.

    Vector<MPI_Request> mpi_requests;
    for (int iproc = 0; iproc < nprocs; ++iproc) {
        if ( ! unique_remote_cols_vv[iproc].empty()) {
            mpi_requests.push_back(MPI_REQUEST_NULL);
            // I need to let other processes know what I need from them.
            BL_MPI_REQUIRE(MPI_Isend(unique_remote_cols_vv[iproc].data(),
                                     int(unique_remote_cols_vv[iproc].size()),
                                     mpi_long, iproc, mpi_tag, mpi_comm,
                                     &(mpi_requests.back())));
            m_recv_from.push_back(iproc);
            m_recv_counts.push_back(int(unique_remote_cols_vv[iproc].size()));
        }
    }

    Vector<Vector<Long>> send_indices(nsends);
    m_total_counts_send = 0;
    for (int isend = 0; isend < nsends; ++isend) {
        MPI_Status mpi_status;
        BL_MPI_REQUIRE(MPI_Probe(MPI_ANY_SOURCE, mpi_tag, mpi_comm, &mpi_status));
        int receiver = mpi_status.MPI_SOURCE;
        int count;
        BL_MPI_REQUIRE(MPI_Get_count(&mpi_status, mpi_long, &count));
        m_send_to.push_back(receiver);
        m_send_counts.push_back(count);
        send_indices[isend].resize(count);
        BL_MPI_REQUIRE(MPI_Recv(send_indices[isend].data(), count, mpi_long,
                                receiver, mpi_tag, mpi_comm, &mpi_status));
        m_total_counts_send += count;
    }

    m_send_indices.resize(m_total_counts_send);
    Gpu::PinnedVector<Long> send_indices_all;
    send_indices_all.reserve(m_total_counts_send);
    for (auto const& vl : send_indices) {
        for (auto x : vl) {
            send_indices_all.push_back(x);
        }
    }
    Gpu::copyAsync(Gpu::hostToDevice, send_indices_all.begin(), send_indices_all.end(),
                   m_send_indices.begin());
    Gpu::streamSynchronize();

    Vector<MPI_Status> mpi_statuses(mpi_requests.size());
    BL_MPI_REQUIRE(MPI_Waitall(int(mpi_requests.size()), mpi_requests.data(),
                               mpi_statuses.data()));

    // Now we convert the remote indices from global to local.
    std::map<Long,Long> gtol;
    for (Long i = 0, N = Long(unique_remote_cols_v.size()); i < N; ++i) {
        gtol[unique_remote_cols_v[i]] = i;
    }
#ifdef AMREX_USE_GPU
    auto& cols = m_remote_cols;
#else
    auto& cols = m_data_remote.col_index;
#endif
    for (auto& c : cols) {
        c = gtol[c];
    }

#ifdef AMREX_USE_GPU
    Gpu::copyAsync(Gpu::hostToDevice, m_remote_cols.begin(), m_remote_cols.end(),
                   m_data_remote.col_index.data());
#endif

    m_comm_prepared = true;
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::pack_buffer (AlgVector<T> const& v)
{
    auto* pdst = m_send_buffer;
    auto* pidx = m_send_indices.data();
    auto const& vv = v.view();
    auto const nsends = Long(m_send_indices.size());
    ParallelFor(nsends, [=] AMREX_GPU_DEVICE (Long i)
    {
        pdst[i] = vv(pidx[i]);
    });
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::unpack_buffer (AlgVector<T>& v)
{
    auto const& csr = m_data_remote;
    if (csr.nnz > 0) {
        T const* AMREX_RESTRICT mat = csr.mat.data();
        auto const* AMREX_RESTRICT col = csr.col_index.data();
        auto const* AMREX_RESTRICT row = csr.row_offset.data();

        auto const* rtol = m_rtol.data();

        auto const* AMREX_RESTRICT px = m_recv_buffer;
        auto      * AMREX_RESTRICT py = v.data();

        auto const nrr = Long(csr.row_offset.size())-1;
        ParallelFor(nrr, [=] AMREX_GPU_DEVICE (Long i)
        {
            T r = 0;
            for (Long j = row[i]; j < row[i+1]; ++j) {
                r += mat[j] * px[col[j]];
            }
            py[rtol[i]] += r;
        });
    }
}

#endif

}

#endif
