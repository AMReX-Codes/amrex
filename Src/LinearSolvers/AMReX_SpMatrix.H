#ifndef AMREX_SP_MATRIX_H_
#define AMREX_SP_MATRIX_H_
#include <AMReX_Config.H>

#include <AMReX_AlgPartition.H>
#include <AMReX_AlgVector.H>
#include <AMReX_CSR.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Gpu.H>
#include <AMReX_Scan.H>

#if defined(AMREX_USE_CUDA)
#  include <cusparse.h>
#elif defined(AMREX_USE_HIP)
#  include <rocsparse/rocsparse.h>
#elif defined(AMREX_USE_SYCL)
#  include <mkl_version.h>
#  include <oneapi/mkl/spblas.hpp>
#endif

#include <fstream>
#include <numeric>
#include <string>
#include <type_traits>
#include <utility>

namespace amrex {

/**
 * \file AMReX_SpMatrix.H
 *
 * Sparse matrix container built on top of CSR partitions (GPU-aware).
 */

//! GPU-ready non-owning CSR data container. The blocks use 32-bit row
//! offsets and column indices: local indices in csr0, indices into
//! `col_map` in csr1.
template <typename T>
struct ParCsr {
    CsrView<T,int> csr0; // diagonal part
    CsrView<T,int> csr1; // off-diagonal part
    Long row_begin = 0; // global row index begin
    Long col_begin = 0; // global col index begin
    Long const* AMREX_RESTRICT row_map = nullptr; // mapping from csr0's row to csr1
    Long const* AMREX_RESTRICT col_map = nullptr; // mapping from csr1's col to global
};

//! Sorted CSR means for each row the column indices are sorted.
struct CsrSorted {
    bool b = true;
    explicit operator bool() const { return b; }
};

//! Valid CSR means all entries are valid. It may be sorted ro unsorted.
struct CsrValid {
    bool b = true;
    explicit operator bool() const { return b; }
};

/**
 * \brief Distributed CSR matrix that manages storage and GPU-friendly partitions.
 */
template <typename T, template<typename> class Allocator = DefaultAllocator>
class SpMatrix
{
public:
    using value_type = T;
    template <class U> using allocator_type = Allocator<U>;
    template <class U> using container_type = PODVector<U,Allocator<U> >;
    using csr_type = CSR<T,container_type>;            // global Long columns (unsplit)
    using local_csr_type = CSR<T,container_type,int>;  // split blocks, 32-bit indices
    using AllocT = Allocator<T>;

    SpMatrix () = default;

    /**
     * \brief Construct a sparse matrix with a fixed number of nonzeros per
     * row.
     *
     * This constructor allocates the required storage and initializes the
     * row offsets. The matrix entries themselves are left
     * uninitialized. Matrix values can be written through the pointer
     * returned by `data()` and the corresponding column indices through the
     * pointer returned by the `columnIndex()`. Alternatively, the matrix
     * may be initialized using `setVal()`. In either case, the CSR data
     * must not contain duplicated entries or invalid entries (i.e.,
     * negative column indices). If the matrix is set manually through
     * pointers and the entries within each row are not sorted by column
     * index, the user must call `sortCSR` after all data have been set.
     *
     * \param partition Parallel partition information.
     * \param nnz_per_row Number of nonzeros per row.
     */
    SpMatrix (AlgPartition partition, int nnz_per_row);

    /**
     * \brief Construct a sparse matrix from a given Partition and CSR
     *
     * The input CSR represents the portion of the matrix local to this
     * process. It must be both sorted and valid. Moreover, it should not
     * have duplicated entries.
     *
     * \param partition Parallel partition information.
     * \param csr CSR represents of the matrix local to this process.
     */
    SpMatrix (AlgPartition partition, csr_type csr);

    SpMatrix (SpMatrix const&) = delete;
    SpMatrix& operator= (SpMatrix const&) = delete;

    SpMatrix (SpMatrix &&) = default;
    SpMatrix& operator= (SpMatrix &&) = default;

    ~SpMatrix () = default;

    /**
     * \brief Allocate storage for a default-constructed matrix with a fixed
     * number of nonzeros per row.
     *
     * This allocates the required storage and initializes the
     * row offsets. The matrix entries themselves are left
     * uninitialized. Matrix values can be written through the pointer
     * returned by `data()` and the corresponding column indices through the
     * pointer returned by the `columnIndex()`. Alternatively, the matrix
     * may be initialized using `setVal()`.
     *
     * \param partition Parallel partition information.
     * \param nnz_per_row Number of nonzeros per row.
     */
    void define (AlgPartition partition, int nnz_per_row);

    /**
     * \brief Define a default-constructed matrix from given CSR arrays.
     *
     * The CSR arrays describe the portion of the matrix local to this
     * process.
     *
     * If `is_valid` indicates that the input CSR may contain invalid
     * entries, this function will filter the input and omit such entries
     * from the resulting matrix. An entry is considered invalid if its
     * column index is negative. In addition, when `is_valid` is false, this
     * function also removes entries whose values are zero. The input CSR is
     * not allowed to contain duplicated valid entries.
     *
     * It is not an error to specify \p is_sorted = false when the input CSR
     * is in fact sorted, nor to specify \p is_valid = false when the input
     * CSR contains no invalid entries. In these cases, the result is
     * correct, but the function may perform unnecessary work and therefore
     * result in a performance penalty. If the input properties are not
     * known, it is safe to set both flags to false.
     *
     * The input data can be freed after this function call. For GPU builds,
     * the data are expected to be in GPU memory.
     *
     * \param partition Parallel partition information.
     * \param mat CSR values array (length `nentries`).
     * \param col_index CSR column indices array (length `nentries`).
     * \param nentries Number of entries.
     * \param row_offset CSR row offsets array (length `num_local_rows + 1`).
     * \param is_sorted Indicates whether the CSR entries are sorted (within each row).
     * \param is_valid Indicates whether the input CSR is guaranteed to be valid
     *                 (i.e., contains no negative column indices).
     */
    void define (AlgPartition partition, T const* mat, Long const* col_index,
                 Long nentries, Long const* row_offset, CsrSorted is_sorted,
                 CsrValid is_valid);

    /**
     * \brief Define a default-constructed matrix from a given CSR.
     *
     * The input CSR represents the portion of the matrix local to this
     * process. It must be valid. Moreover, it should not have duplicated
     * entries.
     *
     * \param partition Parallel partition information.
     * \param csr CSR represents of the matrix local to this process.
     * \param is_sorted Indicates whether the CSR entries are sorted (within each row).
     */
    void define (AlgPartition partition, csr_type csr, CsrSorted is_sorted);

    //! Row partition describing how matrix rows are distributed across ranks.
    [[nodiscard]] AlgPartition const& partition () const { return m_partition; }

    /**
     * \brief Return the column partition used for matrix-vector and
     * matrix-matrix multiplications.
     *
     * This partition corresponds to the row partition of the vector or
     * matrix appearing on the right-hand side of a multiplication.
     */
    [[nodiscard]] AlgPartition const& columnPartition () const { return m_col_partition; }

    //! Number of rows owned by this rank.
    [[nodiscard]] Long numLocalRows () const { return m_row_end - m_row_begin; }
    //! Global row count.
    [[nodiscard]] Long numGlobalRows () const { return m_partition.numGlobalRows(); }
    //! Number of nonzeros stored locally.
    [[nodiscard]] Long numLocalNonZeros () const { return m_nnz; }

    //! Inclusive global index begin on this process.
    [[nodiscard]] Long globalRowBegin () const { return m_row_begin; }
    //! Exclusive global index end on this process.
    [[nodiscard]] Long globalRowEnd () const { return m_row_end; }

    //! Don't use this beyond initial setup.
    [[nodiscard]] T* data () {
        AMREX_ALWAYS_ASSERT(m_split == false);
        return m_csr.mat.data();
    }

    //! Don't use this beyond initial setup.
    [[nodiscard]] Long* columnIndex () {
        AMREX_ALWAYS_ASSERT(m_split == false);
        return m_csr.col_index.data();
    }

    //! Don't use this beyond initial setup.
    [[nodiscard]] Long* rowOffset () {
        AMREX_ALWAYS_ASSERT(m_split == false);
        return m_csr.row_offset.data();
    }

    /*
     * \brief Print the matrix to files.
     *
     * Each process writes its local portion of the matrix to a separate
     * file. This function is provided for debugging purpose. Using it in
     * large-scale runs may overwhelm the file system.
     *
     * File format:
     * - The first line contains three integers describing the local matrix:
     *   the global row begin index, the global row end index, and the
     *   number of local nonzeros.
     * - This is followed by one line per nonzero entry. Each line contains:
     *   the global row index, global column index, and matrix value.
     *
     * \param file Base file name. The full name on process `i` is `{file}.{i}`.
     */
    void printToFile (std::string const& file) const;

    /**
     * \brief Initialize matrix entries using a row-wise functor.
     *
     * The functor must have the signature
     * \code{.cpp}
     * void(Long row, Long* col, T* val)
     * \endcode

     * and is responsible for setting the global column indices and values
     * for the row specified by the global row index. The two pointers point
     * to the beginning of the specified row.
     *
     * This function is intended to be used after
     * \ref define(AlgPartition,int), when the storage has been allocated. It
     * should only be used during initial matrix setup and not for modifying
     * an already assembled matrix.
     *
     * The CSR data must not contain duplicated entries or invalid entries
     * (i.e., negative column indices).
     *
     * \param f Functor `void(Long,Long*,T*)`
     * \param is_sorted Indicates whether the CSR entries are sorted (within each row).
     */
    template <typename F>
    void setVal (F const& f, CsrSorted is_sorted);

    //! Sort CSR by column index withint each row. This is intended to be
    //! used after the user sets the matrix manually or with `setVal`.
    void sortCSR ();

    /**
     * \brief Set the column partition and split the matrix into local and
     * remote blocks with 32-bit column indices.
     *
     * Done on first use by SpMV, SpGEMM and transpose. A matrix keeps its
     * column partition afterwards, and its values can no longer be changed.
     */
    void setColumnPartition (AlgPartition const& col_partition);

    //! Return diagonal elements in a square matrix.
    [[nodiscard]] AlgVector<T,AllocT> const& diagonalVector () const;

    /**
     * \brief Sum the values in each local row and return the result as an AlgVector.
     *
     * Each MPI rank receives the row sums for the rows it owns; if further
     * reductions are needed they should be performed on the returned vector.
     */
    [[nodiscard]] AlgVector<T,AllocT> rowSum () const;

    /**
     * \brief Build GPU-friendly CSR views split into diagonal/off-diagonal blocks.
     *
     * The returned structure aliases the underlying storage; do not outlive this matrix.
     */
    [[nodiscard]] ParCsr<T      >       parcsr ()      ;
    //! Const variant of parcsr().
    [[nodiscard]] ParCsr<T const>       parcsr () const;
    //! Const-qualified alias of parcsr() for convenience.
    [[nodiscard]] ParCsr<T const> const_parcsr () const;

    template <typename U, template<typename> class M, typename N> friend
    void SpMV(AlgVector<U,N>& y, SpMatrix<U,M> const& A, AlgVector<U,N> const& x);

    template <typename U, template<typename> class M> friend
    SpMatrix<U,M> transpose (SpMatrix<U,M> const& A, AlgPartition const& col_partition);

    template <typename U, template<typename> class M> friend
    SpMatrix<U,M> SpGEMM (SpMatrix<U,M> const& A, SpMatrix<U,M> const& B,
                          AlgPartition const& col_partition);

    template <typename U> friend class AMG;

    //! Private helper (exposed for CUDA) that allocates fixed-connectivity matrices with \p nnz_per_row entries.
    void define_doit (int nnz_per_row);

    /**
     * \brief Private helper (exposed for CUDA) that copies/filters CSR arrays into device storage.
     *
     * \tparam I Prefix-scan index type.
     * \param mat        Input values array.
     * \param col_index  Input column indices.
     * \param nentries   Number of candidate entries.
     * \param row_offset Row offsets array.
     */
    template <typename I>
    void define_and_filter_doit (T const* mat, Long const* col_index,
                                 Long nentries, Long const* row_offset);

    /**
     * \brief Prepare halo exchanges for a subsequent SpMV using \p x as the source vector.
     *
     * Must be paired with finishComm_mv(); only available when MPI is enabled.
     */
    void startComm_mv (AlgVector<T,AllocT> const& x);
    //! Finish halo exchanges and accumulate contributions into \p y.
    void finishComm_mv (AlgVector<T,AllocT>& y);

    /**
     * \brief Gather values at the columns of the off-diagonal block.
     *
     * \p x holds local values in the column partition of this matrix
     * (index = global column - column begin). The result is indexed like
     * the column indices of the off-diagonal block, see `parcsr().csr1`.
     * Collective. Without MPI, or when no communication is needed, the
     * result is empty. The column partition must have been set.
     */
    template <typename U>
    Gpu::DeviceVector<U> gatherRemote (U const* x);

    //! Initiate communication required to build the transpose with column partition \p col_partition.
    void startComm_tr (AlgPartition const& col_partition);
    //! Complete transpose communication, writing the assembled matrix into \p AT.
    void finishComm_tr (SpMatrix<T,Allocator>& AT);

public: // NOLINT  Private, but public for CUDA
    //! Split into local and remote blocks with 32-bit column indices.
    void split_csr (AlgPartition const& col_partition);

    /**
     * \brief Define directly in split form from a CSR with compact
     * 32-bit columns: `c < nlocal` is the local column `c`, otherwise
     * `remote_cols[c-nlocal]` is the global column. Rows must be sorted
     * within the local and the remote part.
     */
    void define_split (AlgPartition partition, AlgPartition const& col_partition,
                       local_csr_type&& compact, Long nlocal,
                       Long const* remote_cols, Long nremote);

private:

    void set_num_neighbors ();

    AlgPartition m_partition;
    AlgPartition m_col_partition;
    Long m_row_begin = 0;
    Long m_row_end = 0;
    Long m_col_begin = 0; //! for diagonal part of mat-vec & mat-mat operations
    Long m_col_end = 0;   //! for diagonal part of mat-vec & mat-mat operations
    Long m_nnz = 0; //! total nnz including both diagonal and off-diagonal
    csr_type m_csr;             // unsplit form, empty once split
    local_csr_type m_csr_local; // diagonal block after the split

    mutable AlgVector<T,AllocT> m_diagonal;

    bool m_split = false; // Has the matrix been split into diagonal and off-diagonal parts?

#ifdef AMREX_USE_MPI
    local_csr_type m_csr_remote;

    // The `row_offset` in m_csr_remote does not contain empty rows. The
    // full list of row offsets is available in m_remote_row_offset, which
    // is built on demand by expand_remote_row_offset.
    mutable container_type<int> m_remote_row_offset;

    // It should be noted that m_csr and m_csr_remote may have different
    // number of rows, because some rows may be purely local and they do not
    // appear in m_csr_remote. Thus, we need a mapping from local row index
    // in m_csr_remote to local row index in m_csr. Then we will be able to
    // know its global row index by adding m_row_begin. For example,
    // m_row_begin + m_ri_rtol[i] is the global index for local row i in
    // m_csr_remote.
    container_type<Long> m_ri_rtol; // size: m_csr_remote.nrows()

    // This is local row index mapping from m_csr to m_csr_remote. -1 means
    // the row does not exist in m_csr_remote.
    container_type<Long> m_ri_ltor; // size: numLocalRows()

    // For column index, we also need to be careful with local vs. global,
    // and there two types of locals: local in m_csr and local in
    // m_csr_remote. For m_csr, col_index is the global index if m_split is
    // false, and it becomes the local index if m_split is true. The
    // conversion is global_col_index = local_col_index + m_col_begin.
    //
    // For m_csr_remote, col_index is also local. For a give col_index j,
    // m_remote_cols_v[j] gives us the global index.
    Vector<Long> m_remote_cols_v;
#ifdef AMREX_USE_GPU
    container_type<Long> m_remote_cols_dv;
#endif
    // The size of outer vector is # of procs. The indices are global.
    Vector<Vector<Long>> m_remote_cols_vv;

    // No. of other processes involved in communication.
    //
    // For matrix-vector multiplication, this is the number of processes
    // that require this process' data in the vector. This variable is
    // stored in the matrix and it's set by calling its member function
    // set_num_neighbors(). However, the data here are the vector's data.
    //
    // For matrix transpose, this is the number of processes that will send
    // us their transposed data.
    int m_num_neighbors = -1;

public: // NOLINT  Private functions, but public for cuda

    struct CommMV {
        Vector<int> send_to;
        Vector<int> send_counts;
        Gpu::DeviceVector<Long> send_indices;

        Vector<int> recv_from;
        Vector<int> recv_counts;

        Vector<MPI_Request> send_reqs;
        T* send_buffer = nullptr;
        Long total_counts_send = 0;

        Vector<MPI_Request> recv_reqs;
        T* recv_buffer = nullptr;
        Long total_counts_recv = 0;

        bool prepared = false;
    } m_comm_mv;

    struct CommTR {
        CsrView<T> csrt; // Own the memory inside

        Vector<int> send_to;
        Vector<std::array<int,2>> send_counts;
        Vector<MPI_Request> send_reqs;

        Vector<int> recv_from;
        Vector<std::array<int,2>> recv_counts;
        Vector<MPI_Request> recv_reqs;

        std::array<Long,2> total_counts_recv = {0,0};
        Vector<std::array<Long,4>> recv_buffer_offset;
        // The raw pointers below are owning.
        // xxxxx TODO GPU: Currently we use pinned memory. In the future, we
        // may explore device memory for GPU aware MPI.
        T* recv_buffer_mat = nullptr;
        Long* recv_buffer_col_index = nullptr;
        Long* recv_buffer_row_offset = nullptr;
        Long* recv_buffer_idx_map = nullptr; // local -> global for row of A^T.

    } m_comm_tr;

    //! Builds the remote column list from global column indices and
    //! writes the corresponding 32-bit local indices.
    template <typename VL, typename VI>
    void update_remote_col_index (VL const& gcols, VI& lcols, bool in_device_memory);
    //! Trims the empty rows of the remote block and builds the row maps.
    //! `m_csr_remote.row_offset` must hold the full offsets on entry.
    void trim_remote_rows ();
    void prepare_comm_mv (AlgPartition const& col_partition);
    void pack_buffer_mv (AlgVector<T,AllocT> const& v);
    void unpack_buffer_mv (AlgVector<T,AllocT>& v);
    void comm_tr_recv_wait();
    void comm_tr_clear();
    void unpack_buffer_tr (CommTR const& ctr, AlgPartition const& col_partition);
    void expand_remote_row_offset () const;
    //! Off-diagonal part with the full (untrimmed) row offsets.
    CsrView<T const,int> remote_full_const_view () const;

    //! Rows of another matrix fetched for SpGEMM. Column indices are global.
    struct RemoteRowsMM {
        container_type<Long> row_offset; // size nrows+1
        Long* col_index = nullptr;       // pinned memory
        T* mat = nullptr;                // comms arena
        Long nrows = 0;
        Long nnz = 0;

        RemoteRowsMM () = default;
        ~RemoteRowsMM () { clear(); }
        RemoteRowsMM (RemoteRowsMM const&) = delete;
        RemoteRowsMM& operator= (RemoteRowsMM const&) = delete;
        RemoteRowsMM (RemoteRowsMM&& rhs) noexcept
            : row_offset(std::move(rhs.row_offset)),
              col_index(std::exchange(rhs.col_index, nullptr)),
              mat(std::exchange(rhs.mat, nullptr)),
              nrows(rhs.nrows), nnz(rhs.nnz) {}
        RemoteRowsMM& operator= (RemoteRowsMM&& rhs) noexcept {
            if (this != &rhs) {
                clear();
                row_offset = std::move(rhs.row_offset);
                col_index = std::exchange(rhs.col_index, nullptr);
                mat = std::exchange(rhs.mat, nullptr);
                nrows = rhs.nrows;
                nnz = rhs.nnz;
            }
            return *this;
        }
        void clear () {
            if (col_index || mat) {
                Gpu::streamSynchronize();
                if (col_index) { The_Pinned_Arena()->free(col_index); }
                if (mat) { The_Comms_Arena()->free(mat); }
                col_index = nullptr;
                mat = nullptr;
            }
            row_offset.clear();
            nrows = 0;
            nnz = 0;
        }
    };

    //! Fetch the rows of \p B listed in m_remote_cols_v. This matrix must
    //! already be split with B's row partition. Collective.
    RemoteRowsMM fetch_remote_rows_mm (SpMatrix<T,Allocator> const& B);

#endif
};

namespace detail {

// Communication may be skipped only when all rows and all columns live on
// one common rank.  Two partitions with a single active rank each may still
// name different ranks.
inline bool spmat_comm_is_local (AlgPartition const& row_partition,
                                 AlgPartition const& col_partition)
{
    int const rp = row_partition.singleActiveProc();
    int const cp = col_partition.singleActiveProc();
    return row_partition.numActiveProcs() <= 1
        && col_partition.numActiveProcs() <= 1
        && (rp < 0 || cp < 0 || rp == cp);
}

// Diagonal of a square matrix: entry with column index == row index + offset.
template <typename T, typename I>
void extract_diagonal (T* p, CsrView<T const,I> const& csr, Long offset)
{
    auto const* AMREX_RESTRICT mat = csr.mat;
    auto const* AMREX_RESTRICT col = csr.col_index;
    auto const* AMREX_RESTRICT row = csr.row_offset;
    ParallelForOMP(csr.nrows, [=] AMREX_GPU_DEVICE (Long i)
    {
        T d = 0;
        for (Long j = row[i]; j < row[i+1]; ++j) {
            if (i == Long(col[j]) - offset) {
                d = mat[j];
                break;
            }
        }
        p[i] = d;
    });
}

template <typename T, typename IO, typename II>
void transpose (CsrView<T,IO> const& csrt, CsrView<T const,II> const& csr)
{
    Long nrows = csr.nrows;
    Long ncols = csrt.nrows;
    Long nnz = csr.nnz;

    if (nrows <= 0 || ncols <= 0 || nnz <= 0) {
        auto* p = csrt.row_offset;
        ParallelForOMP(ncols+1, [=] AMREX_GPU_DEVICE (Long i) { p[i] = 0; });
        return;
    }

    // The row indices become the column indices of the transpose.
    if constexpr (!std::is_same_v<IO,Long>) {
        AMREX_ALWAYS_ASSERT(nrows < Long(std::numeric_limits<IO>::max()));
    }

#ifdef AMREX_USE_GPU

#if defined(AMREX_USE_CUDA)

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
        amrex::Abort("SpMatrix transpose: unsupported data type");
    }

    AMREX_ALWAYS_ASSERT(nrows < Long(std::numeric_limits<int>::max()) &&
                        ncols < Long(std::numeric_limits<int>::max()) &&
                        nnz   < Long(std::numeric_limits<int>::max()));

    // cuSPARSE wants 32-bit indices; convert only if needed.
    constexpr bool in32 = std::is_same_v<II,int>;
    constexpr bool out32 = std::is_same_v<IO,int>;
    CsrIndex<int,Gpu::AsyncVector> ci, cit;
    int const* csr_col_index;
    int const* csr_row_offset;
    int* csrt_col_index;
    int* csrt_row_offset;
    if constexpr (in32) {
        csr_col_index = csr.col_index;
        csr_row_offset = csr.row_offset;
    } else {
        ci.copyFrom(csr);
        csr_col_index = ci.col_index.data();
        csr_row_offset = ci.row_offset.data();
    }
    if constexpr (out32) {
        csrt_col_index = csrt.col_index;
        csrt_row_offset = csrt.row_offset;
    } else {
        cit.col_index.resize(csrt.nnz);
        cit.row_offset.resize(csrt.nrows+1);
        csrt_col_index = cit.col_index.data();
        csrt_row_offset = cit.row_offset.data();
    }

    std::size_t buffer_size;
    AMREX_CUSPARSE_SAFE_CALL(
        cusparseCsr2cscEx2_bufferSize(handle, int(nrows), int(ncols), int(nnz),
                                      csr.mat, csr_row_offset, csr_col_index,
                                      csrt.mat, csrt_row_offset, csrt_col_index,
                                      data_type, CUSPARSE_ACTION_NUMERIC,
                                      CUSPARSE_INDEX_BASE_ZERO,
                                      CUSPARSE_CSR2CSC_ALG1,
                                      &buffer_size));

    auto* pbuffer = (void*)The_Async_Arena()->alloc(buffer_size);

    AMREX_CUSPARSE_SAFE_CALL(
        cusparseCsr2cscEx2(handle, int(nrows), int(ncols), int(nnz),
                           csr.mat, csr_row_offset, csr_col_index,
                           csrt.mat, csrt_row_offset, csrt_col_index,
                           data_type, CUSPARSE_ACTION_NUMERIC,
                           CUSPARSE_INDEX_BASE_ZERO,
                           CUSPARSE_CSR2CSC_ALG1,
                           pbuffer));

    if constexpr (!out32) {
        cit.copyTo(csrt);
    }

    // No sync needed: the descriptors are host objects, cusparseDestroy defers
    // the release of GPU resources, and the buffer is freed in stream order.
    AMREX_CUSPARSE_SAFE_CALL(cusparseDestroy(handle));
    The_Async_Arena()->free(pbuffer);

#elif defined(AMREX_USE_HIP)

    rocsparse_handle handle;
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_create_handle(&handle));
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_set_stream(handle, Gpu::gpuStream()));

    constexpr bool same_int = (sizeof(rocsparse_int) == sizeof(II)) &&
                              (sizeof(rocsparse_int) == sizeof(IO));

    rocsparse_int const* csr_col_index;
    rocsparse_int const* csr_row_offset;
    rocsparse_int* csrt_col_index;
    rocsparse_int* csrt_row_offset;
    CsrIndex<rocsparse_int,Gpu::AsyncVector> ci, cit;
    if constexpr (same_int) {
        csr_col_index = reinterpret_cast<rocsparse_int const*>(csr.col_index);
        csr_row_offset = reinterpret_cast<rocsparse_int const*>(csr.row_offset);
        csrt_col_index = reinterpret_cast<rocsparse_int*>(csrt.col_index);
        csrt_row_offset = reinterpret_cast<rocsparse_int*>(csrt.row_offset);
    } else {
        AMREX_ALWAYS_ASSERT(ncols < Long(std::numeric_limits<rocsparse_int>::max()));
        ci.copyFrom(csr);
        cit.col_index.resize(csrt.nnz);
        cit.row_offset.resize(csrt.nrows+1);
        csr_col_index = ci.col_index.data();
        csr_row_offset = ci.row_offset.data();
        csrt_col_index = cit.col_index.data();
        csrt_row_offset = cit.row_offset.data();
    }

    std::size_t buffer_size;
    AMREX_ROCSPARSE_SAFE_CALL(
        rocsparse_csr2csc_buffer_size(handle, rocsparse_int(nrows),
                                      rocsparse_int(ncols), rocsparse_int(nnz),
                                      csr_row_offset, csr_col_index,
                                      rocsparse_action_numeric,
                                      &buffer_size));

    auto* pbuffer = (void*)The_Async_Arena()->alloc(buffer_size);

    if constexpr (std::is_same_v<T,float>) {
        AMREX_ROCSPARSE_SAFE_CALL(
            rocsparse_scsr2csc(handle, rocsparse_int(nrows),
                               rocsparse_int(ncols), rocsparse_int(nnz),
                               csr.mat, csr_row_offset, csr_col_index,
                               csrt.mat, csrt_col_index, csrt_row_offset,
                               rocsparse_action_numeric,
                               rocsparse_index_base_zero,
                               pbuffer));
    } else if constexpr (std::is_same_v<T,double>) {
        AMREX_ROCSPARSE_SAFE_CALL(
            rocsparse_dcsr2csc(handle, rocsparse_int(nrows),
                               rocsparse_int(ncols), rocsparse_int(nnz),
                               csr.mat, csr_row_offset, csr_col_index,
                               csrt.mat, csrt_col_index, csrt_row_offset,
                               rocsparse_action_numeric,
                               rocsparse_index_base_zero,
                               pbuffer));
    } else if constexpr (std::is_same_v<T,GpuComplex<float>>) {
        AMREX_ROCSPARSE_SAFE_CALL(
            rocsparse_ccsr2csc(handle, rocsparse_int(nrows),
                               rocsparse_int(ncols), rocsparse_int(nnz),
                               (rocsparse_float_complex*)csr.mat, csr_row_offset, csr_col_index,
                               (rocsparse_float_complex*)csrt.mat, csrt_col_index, csrt_row_offset,
                               rocsparse_action_numeric,
                               rocsparse_index_base_zero,
                               pbuffer));
    } else if constexpr (std::is_same_v<T,GpuComplex<double>>) {
        AMREX_ROCSPARSE_SAFE_CALL(
            rocsparse_zcsr2csc(handle, rocsparse_int(nrows),
                               rocsparse_int(ncols), rocsparse_int(nnz),
                               (rocsparse_double_complex*)csr.mat, csr_row_offset, csr_col_index,
                               (rocsparse_double_complex*)csrt.mat, csrt_col_index, csrt_row_offset,
                               rocsparse_action_numeric,
                               rocsparse_index_base_zero,
                               pbuffer));
    } else {
        amrex::Abort("SpMatrix transpose: unsupported data type");
    }

    if constexpr (!same_int) {
        cit.copyTo(csrt);
    }

    Gpu::streamSynchronize();
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_destroy_handle(handle));
    The_Async_Arena()->free(pbuffer);

#elif defined(AMREX_USE_SYCL)

    mkl::sparse::matrix_handle_t handle_in{};
    mkl::sparse::matrix_handle_t handle_out{};
    mkl::sparse::init_matrix_handle(&handle_in);
    mkl::sparse::init_matrix_handle(&handle_out);

    // oneMKL needs one index type per handle: Long copies unless already Long.
    Gpu::DeviceVector<Long> lrow_in, lcol_in, lrow_out, lcol_out;
    Long const* prow_in;
    Long const* pcol_in;
    Long* prow_out;
    Long* pcol_out;
    if constexpr (std::is_same_v<II,Long>) {
        prow_in = csr.row_offset;
        pcol_in = csr.col_index;
    } else {
        lrow_in.resize(nrows+1);
        lcol_in.resize(nnz);
        auto* pr = lrow_in.data(); auto const* sr = csr.row_offset;
        auto* pc = lcol_in.data(); auto const* sc = csr.col_index;
        ParallelFor(std::max(nrows+1,nnz), [=] AMREX_GPU_DEVICE (Long i) {
            if (i < nrows+1) { pr[i] = Long(sr[i]); }
            if (i < nnz) { pc[i] = Long(sc[i]); }
        });
        prow_in = lrow_in.data();
        pcol_in = lcol_in.data();
    }
    if constexpr (std::is_same_v<IO,Long>) {
        prow_out = csrt.row_offset;
        pcol_out = csrt.col_index;
    } else {
        lrow_out.resize(ncols+1);
        lcol_out.resize(nnz);
        prow_out = lrow_out.data();
        pcol_out = lcol_out.data();
    }
#if defined(INTEL_MKL_VERSION) && (INTEL_MKL_VERSION < 20250300)
    mkl::sparse::set_csr_data(Gpu::Device::streamQueue(), handle_in, nrows, ncols,
                              mkl::index_base::zero, (Long*)prow_in, (Long*)pcol_in,
                              (T*)csr.mat);
    mkl::sparse::set_csr_data(Gpu::Device::streamQueue(), handle_out, ncols, nrows,
                              mkl::index_base::zero, prow_out, pcol_out, (T*)csrt.mat);
#else
    mkl::sparse::set_csr_data(Gpu::Device::streamQueue(), handle_in, nrows, ncols, nnz,
                              mkl::index_base::zero, (Long*)prow_in, (Long*)pcol_in,
                              (T*)csr.mat);
    mkl::sparse::set_csr_data(Gpu::Device::streamQueue(), handle_out, ncols, nrows, nnz,
                              mkl::index_base::zero, prow_out, pcol_out, (T*)csrt.mat);
#endif

    mkl::sparse::omatcopy(Gpu::Device::streamQueue(), mkl::transpose::trans,
                          handle_in, handle_out);

    mkl::sparse::release_matrix_handle(Gpu::Device::streamQueue(), &handle_in);
    auto ev = mkl::sparse::release_matrix_handle(Gpu::Device::streamQueue(), &handle_out);
    ev.wait();
    if constexpr (!std::is_same_v<IO,Long>) {
        auto* po = csrt.col_index; auto const* so = lcol_out.data();
        auto* pr = csrt.row_offset; auto const* sr = lrow_out.data();
        ParallelFor(std::max(nnz,ncols+1), [=] AMREX_GPU_DEVICE (Long i) {
            if (i < nnz) { po[i] = IO(so[i]); }
            if (i < ncols+1) { pr[i] = IO(sr[i]); }
        });
        Gpu::streamSynchronize();
    }

#endif

    AMREX_GPU_ERROR_CHECK();

#else

    auto* p = csrt.row_offset;

    ParallelForOMP(ncols+1, [=] AMREX_GPU_DEVICE (Long i) { p[i] = 0; });

    // nonzeros per column
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (Long i = 0; i < nnz; ++i) {
        auto col = csr.col_index[i];
#ifdef AMREX_USE_OMP
#pragma omp atomic update
#endif
        ++p[col+1];
    }

    // build row_offset for transposed matrix. Also save a copy.
    Vector<Long> current_pos(ncols+1);
    current_pos[0] = 0;
    for (Long i = 0; i < ncols; ++i) {
        p[i+1] += p[i];
        current_pos[i+1] = p[i+1];
    }

    // The following code is not OMP safe. It's difficult to use OMP and
    // still keep CSR sorted.
    for (Long i = 0; i < nrows; ++i) {
        for (Long idx = csr.row_offset[i]; idx < csr.row_offset[i+1]; ++idx) {
            auto col = csr.col_index[idx];
            Long dest = current_pos[col]++;
            csrt.mat[dest] = csr.mat[idx];
            csrt.col_index[dest] = IO(i);
        }
    }

#endif
}
}

template <typename T, template<typename> class Allocator>
SpMatrix<T,Allocator>::SpMatrix (AlgPartition partition, int nnz_per_row)
    : m_partition(std::move(partition)),
      m_row_begin(m_partition[ParallelDescriptor::MyProc()]),
      m_row_end(m_partition[ParallelDescriptor::MyProc()+1])
{
    define_doit(nnz_per_row);
}

template <typename T, template<typename> class Allocator>
SpMatrix<T,Allocator>::SpMatrix (AlgPartition partition, csr_type csr)
    : m_partition(std::move(partition)),
      m_row_begin(m_partition[ParallelDescriptor::MyProc()]),
      m_row_end(m_partition[ParallelDescriptor::MyProc()+1]),
      m_nnz(csr.nnz),
      m_csr(std::move(csr))
{}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::define (AlgPartition partition, int nnz_per_row)
{
    m_partition = std::move(partition);
    m_row_begin = m_partition[ParallelDescriptor::MyProc()];
    m_row_end = m_partition[ParallelDescriptor::MyProc()+1];
    m_diagonal = AlgVector<T,AllocT>{};
    define_doit(nnz_per_row);
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::define (AlgPartition partition, csr_type csr,
                                    CsrSorted is_sorted)
{
    AMREX_ALWAYS_ASSERT(m_split == false);

    m_partition = std::move(partition);
    m_row_begin = m_partition[ParallelDescriptor::MyProc()];
    m_row_end = m_partition[ParallelDescriptor::MyProc()+1];
    m_nnz = csr.nnz;
    m_csr = std::move(csr);
    m_diagonal = AlgVector<T,AllocT>{};
    if (! is_sorted) { m_csr.sort(); }
}

template <typename T, template<typename> class Allocator>
void
SpMatrix<T,Allocator>::define_doit (int nnz_per_row)
{
    AMREX_ALWAYS_ASSERT(m_split == false);

    nnz_per_row = std::max(nnz_per_row, 0);
    Long nlocalrows = this->numLocalRows();
    m_nnz = nlocalrows*nnz_per_row;
    m_csr.mat.resize(m_nnz);
    m_csr.col_index.resize(m_nnz);
    m_csr.row_offset.resize(nlocalrows+1);
    m_csr.nnz = m_nnz;

    auto* poffset = m_csr.row_offset.data();
    ParallelForOMP(nlocalrows+1, [=] AMREX_GPU_DEVICE (Long lrow) noexcept
    {
        poffset[lrow] = lrow*nnz_per_row;
    });
}

template <typename T, template<typename> class Allocator>
void
SpMatrix<T,Allocator>::define (AlgPartition partition, T const* mat,
                               Long const* col_index, Long nentries,
                               Long const* row_offset, CsrSorted is_sorted,
                               CsrValid is_valid)
{
    AMREX_ALWAYS_ASSERT(m_split == false);

    m_partition = std::move(partition);
    m_row_begin = m_partition[ParallelDescriptor::MyProc()];
    m_row_end = m_partition[ParallelDescriptor::MyProc()+1];
    m_diagonal = AlgVector<T,AllocT>{};

    bool synced = false;

    if (is_valid) {
        m_nnz = nentries;
        Long nlocalrows = this->numLocalRows();
        m_csr.mat.resize(nentries);
        m_csr.col_index.resize(nentries);
        m_csr.row_offset.resize(nlocalrows+1);
        m_csr.nnz = nentries;
        Gpu::copyAsync(Gpu::deviceToDevice, mat, mat+nentries, m_csr.mat.begin());
        Gpu::copyAsync(Gpu::deviceToDevice, col_index, col_index+nentries,
                       m_csr.col_index.begin());
        Gpu::copyAsync(Gpu::deviceToDevice, row_offset, row_offset+nlocalrows+1,
                       m_csr.row_offset.begin());
    } else {
        if (nentries < Long(std::numeric_limits<int>::max())) {
            define_and_filter_doit<int>(mat, col_index, nentries, row_offset);
        } else {
            define_and_filter_doit<Long>(mat, col_index, nentries, row_offset);
        }
        synced = true;
    }

    if (! is_sorted) {
        m_csr.sort();
        synced = true;
    }

    if (! synced) {
        Gpu::streamSynchronize();
    }
}

template <typename T, template<typename> class Allocator>
void
SpMatrix<T,Allocator>::sortCSR ()
{
    AMREX_ALWAYS_ASSERT(m_split == false);
    m_csr.sort();
}

template <typename T, template<typename> class Allocator>
void
SpMatrix<T,Allocator>::setColumnPartition (AlgPartition const& col_partition)
{
    split_csr(col_partition);
}

template <typename T, template<typename> class Allocator>
template <typename I>
void
SpMatrix<T,Allocator>::define_and_filter_doit (T const* mat, Long const* col_index,
                                               Long nentries, Long const* row_offset)
{
    Gpu::DeviceVector<I> psum(nentries);
    auto* ps = psum.data();
    m_nnz = Scan::PrefixSum<I>(I(nentries),
                               [=] AMREX_GPU_DEVICE (I i) -> I {
                                   return col_index[i] >= 0 && mat[i] != 0; },
                               [=] AMREX_GPU_DEVICE (I i, I x) {
                                   ps[i] = x; },
                               Scan::Type::exclusive, Scan::retSum);
    Long nlocalrows = this->numLocalRows();
    m_csr.mat.resize(m_nnz);
    m_csr.col_index.resize(m_nnz);
    m_csr.row_offset.resize(nlocalrows+1);
    m_csr.nnz = m_nnz;
    auto* pmat = m_csr.mat.data();
    auto* pcol = m_csr.col_index.data();
    auto* prow = m_csr.row_offset.data();
    auto actual_nnz = m_nnz;
    ParallelFor(std::max(nentries,nlocalrows+1), [=] AMREX_GPU_DEVICE (Long i)
    {
        if (i < nentries) {
            if (col_index[i] >= 0 && mat[i] != 0) {
                pmat[ps[i]] = mat[i];
                pcol[ps[i]] = col_index[i];
            }
        }
        if (i <= nlocalrows) {
            prow[i] = (i < nlocalrows && row_offset[i] < nentries)
                ? Long(ps[row_offset[i]]) : actual_nnz;
        }
    });
    Gpu::streamSynchronize();
}

template <typename T, template<typename> class Allocator>
void
SpMatrix<T,Allocator>::printToFile (std::string const& file) const
{
    std::ofstream ofs(file+"."+std::to_string(ParallelDescriptor::MyProc()));
    ofs << m_row_begin << " " << m_row_end << " " << m_nnz << "\n";
    Long const nrows = numLocalRows();

    if (!m_split) {
#ifdef AMREX_USE_GPU
        CSR<T,Gpu::PinnedVector> csr;
        amrex::duplicateCSR(Gpu::deviceToHost, csr, m_csr);
        Gpu::streamSynchronize();
#else
        auto const& csr = m_csr;
#endif
        for (Long i = 0; i < nrows; ++i) {
            for (Long j = csr.row_offset[i]; j < csr.row_offset[i+1]; ++j) {
                ofs << i+m_row_begin << " " << csr.col_index[j] << " " << csr.mat[j] << "\n";
            }
        }
        return;
    }

#ifdef AMREX_USE_GPU
    CSR<T,Gpu::PinnedVector,int> csr;
    amrex::duplicateCSR(Gpu::deviceToHost, csr, m_csr_local);
#  ifdef AMREX_USE_MPI
    CSR<T,Gpu::PinnedVector,int> csr_r;
    amrex::duplicateCSR(Gpu::deviceToHost, csr_r, m_csr_remote);
    Gpu::PinnedVector<Long> ri_ltor(m_ri_ltor.size());
    Gpu::copyAsync(Gpu::deviceToHost, m_ri_ltor.begin(), m_ri_ltor.end(), ri_ltor.begin());
#  endif
    Gpu::streamSynchronize();
#else
    auto const& csr = m_csr_local;
#  ifdef AMREX_USE_MPI
    auto const& csr_r = m_csr_remote;
    auto const& ri_ltor = m_ri_ltor;
#  endif
#endif

    for (Long i = 0; i < nrows; ++i) {
        for (Long j = csr.row_offset[i]; j < csr.row_offset[i+1]; ++j) {
            ofs << i+m_row_begin << " " << Long(csr.col_index[j])+m_col_begin
                << " " << csr.mat[j] << "\n";
        }
#ifdef AMREX_USE_MPI
        if (i < Long(ri_ltor.size()) && ri_ltor[i] >= 0) {
            Long ii = ri_ltor[i];
            for (Long j = csr_r.row_offset[ii]; j < csr_r.row_offset[ii+1]; ++j) {
                ofs << i+m_row_begin << " " << m_remote_cols_v[csr_r.col_index[j]]
                    << " " << csr_r.mat[j] << "\n";
            }
        }
#endif
    }
}

template <typename T, template<typename> class Allocator>
template <typename F>
void SpMatrix<T,Allocator>::setVal (F const& f, CsrSorted is_sorted)
{
    // xxxxx TODO: We can try to optimize this later by using shared memory.

    AMREX_ALWAYS_ASSERT(m_split == false);
    m_diagonal = AlgVector<T,AllocT>{};

    Long nlocalrows = this->numLocalRows();
    Long rowbegin = this->globalRowBegin();
    auto* pmat = m_csr.mat.data();
    auto* pcolindex = m_csr.col_index.data();
    auto* prowoffset = m_csr.row_offset.data();
    ParallelForOMP(nlocalrows, [=] AMREX_GPU_DEVICE (int lrow) noexcept
    {
        f(rowbegin+lrow, pcolindex+prowoffset[lrow], pmat+prowoffset[lrow]);
    });

    if (! is_sorted) { m_csr.sort(); }
}

template <typename T, template<typename> class Allocator>
AlgVector<T,Allocator<T>> const& SpMatrix<T,Allocator>::diagonalVector () const
{
    if (m_diagonal.empty()) {
        // Square matrix: global columns before the split, local block after.
        if (m_split) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_col_partition == this->partition(),
                                             "SpMatrix::diagonalVector: column partition differs from row partition");
        }
        m_diagonal.define(this->partition());
        if (m_split) {
            detail::extract_diagonal(m_diagonal.data(), m_csr_local.const_view(), Long(0));
        } else {
            detail::extract_diagonal(m_diagonal.data(), m_csr.const_view(), m_row_begin);
        }
    }
    return m_diagonal;
}

template <typename T, template<typename> class Allocator>
AlgVector<T,Allocator<T>> SpMatrix<T,Allocator>::rowSum () const
{
    AlgVector<T,Allocator<T>> r(this->partition());
    auto* p = r.data();
    if (!m_split) {
        auto const c = m_csr.const_view();
        ParallelForOMP(r.numLocalRows(), [=] AMREX_GPU_DEVICE (Long i)
        {
            T s = 0;
            for (auto idx = c.row_offset[i]; idx < c.row_offset[i+1]; ++idx) {
                s += c.mat[idx];
            }
            p[i] = s;
        });
        return r;
    }
    auto const& a = this->const_parcsr();
    ParallelForOMP(r.numLocalRows(), [=] AMREX_GPU_DEVICE (Long i)
    {
        T s = 0;
        for (auto idx = a.csr0.row_offset[i];
                  idx < a.csr0.row_offset[i+1]; ++idx) {
            s += a.csr0.mat[idx];
        }
        if (a.csr1.nnz > 0 && a.row_map[i] >= 0) {
            auto ii = a.row_map[i];
            for (auto idx = a.csr1.row_offset[ii];
                      idx < a.csr1.row_offset[ii+1]; ++idx) {
                s += a.csr1.mat[idx];
            }
        }
        p[i] = s;
    });
    return r;
}

template <typename T, template<typename> class Allocator>
ParCsr<T> SpMatrix<T,Allocator>::parcsr ()
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_split, "SpMatrix::parcsr: column partition not set");
    return ParCsr<T>{m_csr_local.view(),
#ifdef AMREX_USE_MPI
                     m_csr_remote.view(),
#else
                     CsrView<T,int>{},
#endif
                     m_row_begin,
                     m_col_begin,
#ifdef AMREX_USE_MPI
                     m_ri_ltor.data(),
#  ifdef AMREX_USE_GPU
                     m_remote_cols_dv.data()
#  else
                     m_remote_cols_v.data()
#  endif
#else
                     nullptr, nullptr
#endif
                    };
}

template <typename T, template<typename> class Allocator>
ParCsr<T const> SpMatrix<T,Allocator>::const_parcsr () const
{
    using U = T const;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_split, "SpMatrix::parcsr: column partition not set");
    return ParCsr<U>{m_csr_local.const_view(),
#ifdef AMREX_USE_MPI
                     m_csr_remote.const_view(),
#else
                     CsrView<U,int>{},
#endif
                     m_row_begin,
                     m_col_begin,
#ifdef AMREX_USE_MPI
                     m_ri_ltor.data(),
#  ifdef AMREX_USE_GPU
                     m_remote_cols_dv.data()
#  else
                     m_remote_cols_v.data()
#  endif
#else
                     nullptr, nullptr
#endif
                    };
}

template <typename T, template<typename> class Allocator>
ParCsr<T const> SpMatrix<T,Allocator>::parcsr () const
{
    return this->const_parcsr();
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::startComm_mv (AlgVector<T,AllocT> const& x)
{
    BL_PROFILE("SpMatrix::startComm_mv");
#ifndef AMREX_USE_MPI
    amrex::ignore_unused(this, x);
#else
    if (detail::spmat_comm_is_local(this->partition(), x.partition())) { return; }

    this->prepare_comm_mv(x.partition());

    auto const mpi_tag = ParallelDescriptor::SeqNum();
    auto const mpi_t_type = ParallelDescriptor::Mpi_typemap<T>::type();
    auto const mpi_comm = ParallelContext::CommunicatorSub();

    auto const nrecvs = int(m_comm_mv.recv_from.size());
    if (nrecvs > 0) {
        m_comm_mv.recv_buffer = (T*)The_Comms_Arena()->alloc(sizeof(T)*m_comm_mv.total_counts_recv);
        m_comm_mv.recv_reqs.resize(nrecvs, MPI_REQUEST_NULL);
        auto* p_recv = m_comm_mv.recv_buffer;
        for (int irecv = 0; irecv < nrecvs; ++irecv) {
            BL_MPI_REQUIRE(MPI_Irecv(p_recv,
                                     m_comm_mv.recv_counts[irecv], mpi_t_type,
                                     m_comm_mv.recv_from[irecv], mpi_tag, mpi_comm,
                                     &(m_comm_mv.recv_reqs[irecv])));
            p_recv += m_comm_mv.recv_counts[irecv];
        }
        AMREX_ASSERT(p_recv == m_comm_mv.recv_buffer + m_comm_mv.total_counts_recv);
    }

    auto const nsends = int(m_comm_mv.send_to.size());
    if (nsends > 0) {
        m_comm_mv.send_buffer = (T*)The_Comms_Arena()->alloc(sizeof(T)*m_comm_mv.total_counts_send);

        pack_buffer_mv(x);
        Gpu::streamSynchronize();

        m_comm_mv.send_reqs.resize(nsends, MPI_REQUEST_NULL);
        auto* p_send = m_comm_mv.send_buffer;
        for (int isend = 0; isend < nsends; ++isend) {
            auto count = m_comm_mv.send_counts[isend];
            BL_MPI_REQUIRE(MPI_Isend(p_send, count, mpi_t_type, m_comm_mv.send_to[isend],
                                     mpi_tag, mpi_comm, &(m_comm_mv.send_reqs[isend])));
            p_send += count;
        }
        AMREX_ASSERT(p_send == m_comm_mv.send_buffer + m_comm_mv.total_counts_send);
    }
#endif
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::finishComm_mv (AlgVector<T,AllocT>& y)
{
    BL_PROFILE("SpMatrix::finishComm_mv");
#ifndef AMREX_USE_MPI
    amrex::ignore_unused(y);
#else
    if (detail::spmat_comm_is_local(this->partition(), m_col_partition)) { return; }

    if ( ! m_comm_mv.recv_reqs.empty()) {
        Vector<MPI_Status> mpi_statuses(m_comm_mv.recv_reqs.size());
        BL_MPI_REQUIRE(MPI_Waitall(int(m_comm_mv.recv_reqs.size()),
                                   m_comm_mv.recv_reqs.data(),
                                   mpi_statuses.data()));
    }

    unpack_buffer_mv(y);

    if ( ! m_comm_mv.send_reqs.empty()) {
        Vector<MPI_Status> mpi_statuses(m_comm_mv.send_reqs.size());
        BL_MPI_REQUIRE(MPI_Waitall(int(m_comm_mv.send_reqs.size()),
                                   m_comm_mv.send_reqs.data(),
                                   mpi_statuses.data()));
    }

    Gpu::streamSynchronize();
    The_Comms_Arena()->free(m_comm_mv.send_buffer);
    The_Comms_Arena()->free(m_comm_mv.recv_buffer);
    m_comm_mv.send_reqs.clear();
    m_comm_mv.recv_reqs.clear();
#endif
}

template <typename T, template<typename> class Allocator>
template <typename U>
Gpu::DeviceVector<U> SpMatrix<T,Allocator>::gatherRemote (U const* x)
{
    BL_PROFILE("SpMatrix::gatherRemote");
    Gpu::DeviceVector<U> r;
#ifndef AMREX_USE_MPI
    amrex::ignore_unused(this, x);
#else
    AMREX_ALWAYS_ASSERT(m_split);
    if (detail::spmat_comm_is_local(this->partition(), m_col_partition)) { return r; }

    this->prepare_comm_mv(m_col_partition);

    auto const mpi_tag = ParallelDescriptor::SeqNum();
    auto const mpi_type = ParallelDescriptor::Mpi_typemap<U>::type();
    auto const mpi_comm = ParallelContext::CommunicatorSub();

    auto const nrecvs = int(m_comm_mv.recv_from.size());
    U* recv_buffer = nullptr;
    Vector<MPI_Request> recv_reqs(nrecvs, MPI_REQUEST_NULL);
    if (nrecvs > 0) {
        recv_buffer = (U*)The_Comms_Arena()->alloc(sizeof(U)*m_comm_mv.total_counts_recv);
        auto* p_recv = recv_buffer;
        for (int irecv = 0; irecv < nrecvs; ++irecv) {
            BL_MPI_REQUIRE(MPI_Irecv(p_recv, m_comm_mv.recv_counts[irecv], mpi_type,
                                     m_comm_mv.recv_from[irecv], mpi_tag, mpi_comm,
                                     &recv_reqs[irecv]));
            p_recv += m_comm_mv.recv_counts[irecv];
        }
    }

    auto const nsends = int(m_comm_mv.send_to.size());
    U* send_buffer = nullptr;
    Vector<MPI_Request> send_reqs(nsends, MPI_REQUEST_NULL);
    if (nsends > 0) {
        send_buffer = (U*)The_Comms_Arena()->alloc(sizeof(U)*m_comm_mv.total_counts_send);
        auto* AMREX_RESTRICT pdst = send_buffer;
        auto const* AMREX_RESTRICT pidx = m_comm_mv.send_indices.data();
        auto const col_begin = m_col_begin;
        ParallelForOMP(m_comm_mv.total_counts_send, [=] AMREX_GPU_DEVICE (Long i)
        {
            pdst[i] = x[pidx[i]-col_begin];
        });
        Gpu::streamSynchronize();

        auto* p_send = send_buffer;
        for (int isend = 0; isend < nsends; ++isend) {
            auto count = m_comm_mv.send_counts[isend];
            BL_MPI_REQUIRE(MPI_Isend(p_send, count, mpi_type, m_comm_mv.send_to[isend],
                                     mpi_tag, mpi_comm, &send_reqs[isend]));
            p_send += count;
        }
    }

    r.resize(m_comm_mv.total_counts_recv);
    if (nrecvs > 0) {
        Vector<MPI_Status> mpi_statuses(nrecvs);
        BL_MPI_REQUIRE(MPI_Waitall(nrecvs, recv_reqs.data(), mpi_statuses.data()));
        // The buffer may be pinned host memory, so copy with a kernel.
        auto* AMREX_RESTRICT pr = r.data();
        auto const* AMREX_RESTRICT pb = recv_buffer;
        ParallelForOMP(m_comm_mv.total_counts_recv, [=] AMREX_GPU_DEVICE (Long i)
        {
            pr[i] = pb[i];
        });
    }
    if (nsends > 0) {
        Vector<MPI_Status> mpi_statuses(nsends);
        BL_MPI_REQUIRE(MPI_Waitall(nsends, send_reqs.data(), mpi_statuses.data()));
    }

    Gpu::streamSynchronize();
    The_Comms_Arena()->free(send_buffer);
    The_Comms_Arena()->free(recv_buffer);
#endif
    return r;
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::startComm_tr (AlgPartition const& col_partition)
{
    BL_PROFILE("SpMatrix::startComm_tr");
#ifdef AMREX_USE_MPI
    if (detail::spmat_comm_is_local(this->partition(), col_partition)) { return; }

    this->split_csr(col_partition);

    int const nprocs = ParallelContext::NProcsSub();
    auto const mpi_tag = ParallelDescriptor::SeqNum();
    auto const mpi_long = ParallelDescriptor::Mpi_typemap<Long>::type();
    auto const mpi_t = ParallelDescriptor::Mpi_typemap<T>::type();
    auto const mpi_comm = ParallelContext::CommunicatorSub();

    // transpose the off-diagonal part
    if (m_csr_remote.nnz > 0) {
        m_comm_tr.csrt.nnz = m_csr_remote.nnz;
        m_comm_tr.csrt.nrows = m_remote_cols_v.size();
        m_comm_tr.csrt.mat = (T*)The_Pinned_Arena()->alloc
            (sizeof(T)*m_comm_tr.csrt.nnz);
        m_comm_tr.csrt.col_index = (Long*)The_Pinned_Arena()->alloc
            (sizeof(Long)*m_comm_tr.csrt.nnz);
        m_comm_tr.csrt.row_offset = (Long*)The_Pinned_Arena()->alloc
            (sizeof(Long)*(m_comm_tr.csrt.nrows+1));
#ifdef AMREX_USE_GPU
        csr_type csr_comm;
        csr_comm.resize(m_comm_tr.csrt.nrows, m_comm_tr.csrt.nnz);
        auto const& csrv_comm = csr_comm.view();
#else
        auto const& csrv_comm = m_comm_tr.csrt;
#endif
        detail::transpose(csrv_comm, m_csr_remote.const_view());
        auto row_begin = m_row_begin;
        auto ri_rtol = m_ri_rtol.data();
        auto* col_index = csrv_comm.col_index;
        ParallelForOMP(csrv_comm.nnz, [=] AMREX_GPU_DEVICE (Long idx)
        {
            auto gjt =ri_rtol[col_index[idx]] + row_begin;
            col_index[idx] = gjt; // global index
        });
#ifdef AMREX_USE_GPU
        Gpu::copyAsync(Gpu::deviceToHost,
                       csrv_comm.     mat,
                       csrv_comm.     mat + csrv_comm.nnz,
                       m_comm_tr.csrt.mat);
        Gpu::copyAsync(Gpu::deviceToHost,
                       csrv_comm.     col_index,
                       csrv_comm.     col_index + csrv_comm.nnz,
                       m_comm_tr.csrt.col_index);
        Gpu::copyAsync(Gpu::deviceToHost,
                       csrv_comm.     row_offset,
                       csrv_comm.     row_offset + csrv_comm.nrows+1,
                       m_comm_tr.csrt.row_offset);
        Gpu::streamSynchronize();
#endif
    }

    if (m_num_neighbors < 0) { set_num_neighbors(); }

    // As a sender, I need to let other processes know that how many
    // elements I will send them.

    Vector<MPI_Request> mpi_requests;
    mpi_requests.reserve(nprocs);
    if (m_csr_remote.nnz > 0) {
        Long it = 0;
        for (int iproc = 0; iproc < nprocs; ++iproc) {
            Long n = 0;
            for (Long i = 0; i < Long(m_remote_cols_vv[iproc].size()); ++i) {
                n += m_comm_tr.csrt.row_offset[it+1] - m_comm_tr.csrt.row_offset[it];
                ++it;
            }
            if (n > 0) {
                mpi_requests.push_back(MPI_REQUEST_NULL);
                AMREX_ALWAYS_ASSERT(n < std::numeric_limits<int>::max());
                std::array<int,2> nn{int(n), int(m_remote_cols_vv[iproc].size())};
                BL_MPI_REQUIRE(MPI_Isend(nn.data(), 2, MPI_INT, iproc, mpi_tag,
                                         mpi_comm, &(mpi_requests.back())));
                m_comm_tr.send_to.push_back(iproc);
                m_comm_tr.send_counts.push_back(nn);
            }
        }
    }

    // As a receiver, m_num_neighbors is the number of processes from which
    // I will receive data.

    for (int irecv = 0; irecv < m_num_neighbors; ++irecv) {
        MPI_Status mpi_status;
        BL_MPI_REQUIRE(MPI_Probe(MPI_ANY_SOURCE, mpi_tag, mpi_comm, &mpi_status));
        int sender = mpi_status.MPI_SOURCE;
        std::array<int,2> nn;
        BL_MPI_REQUIRE(MPI_Recv(nn.data(), 2, MPI_INT, sender, mpi_tag,
                                mpi_comm, &mpi_status));
        m_comm_tr.recv_from.push_back(sender);
        m_comm_tr.recv_counts.push_back(nn);
        m_comm_tr.total_counts_recv[0] += nn[0];
        m_comm_tr.total_counts_recv[1] += nn[1];
    }

    if (! mpi_requests.empty()) {
        Vector<MPI_Status> mpi_statuses(mpi_requests.size());
        BL_MPI_REQUIRE(MPI_Waitall(int(mpi_requests.size()), mpi_requests.data(),
                                   mpi_statuses.data()));
    }

    auto const mpi_tag_m = ParallelDescriptor::SeqNum();
    auto const mpi_tag_c = ParallelDescriptor::SeqNum();
    auto const mpi_tag_r = ParallelDescriptor::SeqNum();
    auto const mpi_tag_p = ParallelDescriptor::SeqNum();

    // We need to send m_comm_tr.csrt.mat, col_index & row_offset. We also
    // need to send m_remote_cols_vv, which maps row index (in transposed
    // matrix) form local to global.

    auto const nrecvs = int(m_comm_tr.recv_from.size());
    if (nrecvs > 0) {
        m_comm_tr.recv_buffer_mat = (T*) The_Pinned_Arena()->alloc
            (sizeof(T) * m_comm_tr.total_counts_recv[0]);
        m_comm_tr.recv_buffer_col_index = (Long*) The_Pinned_Arena()->alloc
            (sizeof(Long) * m_comm_tr.total_counts_recv[0]);
        m_comm_tr.recv_buffer_row_offset = (Long*) The_Pinned_Arena()->alloc
            (sizeof(Long) * (m_comm_tr.total_counts_recv[1]+nrecvs));
        m_comm_tr.recv_buffer_idx_map = (Long*) The_Pinned_Arena()->alloc
            (sizeof(Long) * m_comm_tr.total_counts_recv[1]);
        m_comm_tr.recv_buffer_offset.push_back({0,0,0,0});
        m_comm_tr.recv_reqs.resize(4*nrecvs, MPI_REQUEST_NULL);
        for (int irecv = 0; irecv < nrecvs; ++irecv) {
            auto [os0, os1, os2, os3] = m_comm_tr.recv_buffer_offset.back();
            auto [n0, n1] = m_comm_tr.recv_counts[irecv];
            auto recv_from_rank = m_comm_tr.recv_from[irecv];
            BL_MPI_REQUIRE(MPI_Irecv(m_comm_tr.recv_buffer_mat + os0,
                                     n0,
                                     mpi_t,
                                     recv_from_rank,
                                     mpi_tag_m,
                                     mpi_comm,
                                     &(m_comm_tr.recv_reqs[irecv*4])));
            BL_MPI_REQUIRE(MPI_Irecv(m_comm_tr.recv_buffer_col_index + os1,
                                     n0,
                                     mpi_long,
                                     recv_from_rank,
                                     mpi_tag_c,
                                     mpi_comm,
                                     &(m_comm_tr.recv_reqs[irecv*4+1])));
            BL_MPI_REQUIRE(MPI_Irecv(m_comm_tr.recv_buffer_row_offset + os2,
                                     n1+1,
                                     mpi_long,
                                     recv_from_rank,
                                     mpi_tag_r,
                                     mpi_comm,
                                     &(m_comm_tr.recv_reqs[irecv*4+2])));
            BL_MPI_REQUIRE(MPI_Irecv(m_comm_tr.recv_buffer_idx_map + os3,
                                     n1,
                                     mpi_long,
                                     recv_from_rank,
                                     mpi_tag_p,
                                     mpi_comm,
                                     &(m_comm_tr.recv_reqs[irecv*4+3])));
            m_comm_tr.recv_buffer_offset.push_back({os0 + n0,
                                                    os1 + n0,
                                                    os2 + n1+1,
                                                    os3 + n1});
        }
    }

    auto const nsends = int(m_comm_tr.send_to.size());
    if (nsends > 0) {
        m_comm_tr.send_reqs.resize(4*nsends, MPI_REQUEST_NULL);
        Long os0 = 0, os1 = 0;
        for (int isend = 0; isend < nsends; ++isend) {
            auto [n0, n1] = m_comm_tr.send_counts[isend];
            auto send_to_rank = m_comm_tr.send_to[isend];
            BL_MPI_REQUIRE(MPI_Isend(m_comm_tr.csrt.mat + os0,
                                     n0,
                                     mpi_t,
                                     send_to_rank,
                                     mpi_tag_m,
                                     mpi_comm,
                                     &(m_comm_tr.send_reqs[isend*4])));
            BL_MPI_REQUIRE(MPI_Isend(m_comm_tr.csrt.col_index + os0,
                                     n0,
                                     mpi_long,
                                     send_to_rank,
                                     mpi_tag_c,
                                     mpi_comm,
                                     &(m_comm_tr.send_reqs[isend*4+1])));
            BL_MPI_REQUIRE(MPI_Isend(m_comm_tr.csrt.row_offset + os1,
                                     n1+1,
                                     mpi_long,
                                     send_to_rank,
                                     mpi_tag_r,
                                     mpi_comm,
                                     &(m_comm_tr.send_reqs[isend*4+2])));
            BL_MPI_REQUIRE(MPI_Isend(m_remote_cols_vv[send_to_rank].data(),
                                     n1,
                                     mpi_long,
                                     send_to_rank,
                                     mpi_tag_p,
                                     mpi_comm,
                                     &(m_comm_tr.send_reqs[isend*4+3])));
            os0 += n0;
            os1 += n1;
        }
    }
#else
    amrex::ignore_unused(col_partition);
#endif
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::finishComm_tr (SpMatrix<T,Allocator>& AT)
{
    BL_PROFILE("SpMatrix::finishComm_tr");
#ifdef AMREX_USE_MPI
    if (detail::spmat_comm_is_local(this->partition(), AT.partition())) { return; }

    this->comm_tr_recv_wait();

    AT.unpack_buffer_tr(m_comm_tr, this->m_partition);

    this->comm_tr_clear();
#else
    amrex::ignore_unused(AT);
#endif
}

#ifdef AMREX_USE_MPI

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::comm_tr_recv_wait ()
{
    if (! m_comm_tr.recv_reqs.empty()) {
        Vector<MPI_Status> mpi_statuses(m_comm_tr.recv_reqs.size());
        BL_MPI_REQUIRE(MPI_Waitall(int(m_comm_tr.recv_reqs.size()),
                                   m_comm_tr.recv_reqs.data(),
                                   mpi_statuses.data()));
    }
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::comm_tr_clear ()
{
    if (! m_comm_tr.send_reqs.empty()) {
        Vector<MPI_Status> mpi_statuses(m_comm_tr.send_reqs.size());
        BL_MPI_REQUIRE(MPI_Waitall(int(m_comm_tr.send_reqs.size()),
                                   m_comm_tr.send_reqs.data(),
                                   mpi_statuses.data()));
    }

    if (m_comm_tr.csrt.nnz > 0) {
        The_Pinned_Arena()->free(m_comm_tr.csrt.mat);
        The_Pinned_Arena()->free(m_comm_tr.csrt.col_index);
        The_Pinned_Arena()->free(m_comm_tr.csrt.row_offset);
    }
    if (m_comm_tr.recv_buffer_mat) {
        The_Pinned_Arena()->free(m_comm_tr.recv_buffer_mat);
        The_Pinned_Arena()->free(m_comm_tr.recv_buffer_col_index);
        The_Pinned_Arena()->free(m_comm_tr.recv_buffer_row_offset);
        The_Pinned_Arena()->free(m_comm_tr.recv_buffer_idx_map);
    }
    m_comm_tr = CommTR{};
}

#endif

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::define_split (AlgPartition partition,
                                          AlgPartition const& col_partition,
                                          local_csr_type&& compact, Long nlocal,
                                          Long const* remote_cols, Long nremote)
{
    BL_PROFILE("SpMatrix::define_split");
    AMREX_ALWAYS_ASSERT(m_split == false && m_col_partition.empty());
    amrex::ignore_unused(nremote);

    m_partition = std::move(partition);
    m_row_begin = m_partition[ParallelDescriptor::MyProc()];
    m_row_end = m_partition[ParallelDescriptor::MyProc()+1];
    AMREX_ASSERT(compact.nrows() == this->numLocalRows());
    AMREX_ASSERT(nremote >= 0 && (nremote == 0 || remote_cols != nullptr));
    m_diagonal = AlgVector<T,AllocT>{};
    m_col_partition = col_partition;
    m_col_begin = col_partition[ParallelDescriptor::MyProc()];
    m_col_end   = col_partition[ParallelDescriptor::MyProc()+1];
    AMREX_ALWAYS_ASSERT(m_col_end - m_col_begin == nlocal);
    m_nnz = compact.nnz;
    m_csr = csr_type{};

    Long const nlocalrows = this->numLocalRows();
    Long const nnz = compact.nnz;

    // Entries with a compact column below nlocal are local; within a row
    // they come first.
    Gpu::DeviceVector<Long> pfsum(nnz);
    auto* p_pfsum = pfsum.data();
    auto const* pcol = compact.col_index.data();
    AMREX_ALWAYS_ASSERT(nlocal < Long(std::numeric_limits<int>::max()));
    int const nl = int(nlocal);
    Long const local_nnz = Scan::PrefixSum<Long>(nnz,
        [=] AMREX_GPU_DEVICE (Long i) -> Long { return pcol[i] < nl; },
        [=] AMREX_GPU_DEVICE (Long i, Long const& x) { p_pfsum[i] = x; },
        Scan::Type::exclusive, Scan::retSum);
    Long const remote_nnz = nnz - local_nnz;

#ifndef AMREX_USE_MPI
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(remote_nnz == 0,
        "SpMatrix::define_split: remote columns without MPI");
#endif

    if (remote_nnz == 0) {
        m_csr_local = std::move(compact);
    } else {
        m_csr_local.resize(nlocalrows, local_nnz);
        container_type<Long> remote_gcols(remote_nnz);
        container_type<T> remote_mat(remote_nnz);
        container_type<int> remote_row_offset(nlocalrows+1);
        auto const* pmat = compact.mat.data();
        auto* pmat_l = m_csr_local.mat.data();
        auto* pcol_l = m_csr_local.col_index.data();
        auto* pmat_r = remote_mat.data();
        auto* pcol_r = remote_gcols.data();
        ParallelForOMP(nnz, [=] AMREX_GPU_DEVICE (Long i)
        {
            auto ps = p_pfsum[i];
            if (pcol[i] < nl) {
                pmat_l[ps] = pmat[i];
                pcol_l[ps] = pcol[i];
            } else {
                pmat_r[i-ps] = pmat[i];
                pcol_r[i-ps] = remote_cols[pcol[i] - nl];
            }
        });
        auto const noffset = nlocalrows+1;
        auto const* pro = compact.row_offset.data();
        auto* pro_l = m_csr_local.row_offset.data();
        auto* pro_r = remote_row_offset.data();
        ParallelForOMP(noffset, [=] AMREX_GPU_DEVICE (Long i)
        {
            Long ro_l = (i < noffset-1)
                ? ((pro[i] < nnz) ? p_pfsum[pro[i]] : local_nnz) : local_nnz;
            pro_l[i] = int(ro_l);
            pro_r[i] = int(pro[i] - ro_l);
        });
        Gpu::streamSynchronize();
        compact = local_csr_type{};
#ifdef AMREX_USE_MPI
        m_csr_remote.mat = std::move(remote_mat);
        m_csr_remote.col_index.resize(remote_nnz);
        m_csr_remote.row_offset = std::move(remote_row_offset);
        m_csr_remote.nnz = remote_nnz;
        trim_remote_rows();
        update_remote_col_index(remote_gcols, m_csr_remote.col_index, true);
#endif
    }
#ifdef AMREX_USE_MPI
    if (remote_nnz == 0) {
        // Keep the remote bookkeeping consistent with an empty block.
        update_remote_col_index(container_type<Long>{}, m_csr_remote.col_index, true);
    }
#endif

    m_split = true;
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::split_csr (AlgPartition const& col_partition)
{
    BL_PROFILE("SpMatrix::split_csr");
    if (m_split) {
        AMREX_ALWAYS_ASSERT
            (m_col_begin == col_partition[ParallelDescriptor::MyProc()] &&
             m_col_end   == col_partition[ParallelDescriptor::MyProc()+1]);
        return;
    }

    AMREX_ALWAYS_ASSERT(m_col_partition.empty());

    m_col_partition = col_partition;
    m_col_begin = col_partition[ParallelDescriptor::MyProc()];
    m_col_end   = col_partition[ParallelDescriptor::MyProc()+1];
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE
        (m_col_end - m_col_begin < Long(std::numeric_limits<int>::max()),
         "SpMatrix: the local column block is too large for 32-bit indices");

    // This function needs to be safe when nnz is zero.

    // Split into a diagonal block with 32-bit local column indices and,
    // with MPI, an off-diagonal block with 32-bit indices into the sorted
    // list of remote columns.

    Long const nlocalrows = this->numLocalRows();
    Long local_nnz;
    Gpu::DeviceVector<Long> pfsum(m_nnz);
    auto* p_pfsum = pfsum.data();
    auto col_begin = m_col_begin;
    auto col_end = m_col_end;
    auto const* pcol = m_csr.col_index.data();
    if (m_csr.nnz < Long(std::numeric_limits<int>::max())) {
        local_nnz = Scan::PrefixSum<int>(int(m_nnz),
                                         [=] AMREX_GPU_DEVICE (int i) -> int {
                                             return (pcol[i] >= col_begin &&
                                                     pcol[i] <  col_end); },
                                         [=] AMREX_GPU_DEVICE (int i, int const& x) {
                                             p_pfsum[i] = x; },
                                         Scan::Type::exclusive, Scan::retSum);
    } else {
        local_nnz = Scan::PrefixSum<Long>(m_nnz,
                                          [=] AMREX_GPU_DEVICE (Long i) -> Long {
                                              return (pcol[i] >= col_begin &&
                                                      pcol[i] <  col_end); },
                                          [=] AMREX_GPU_DEVICE (Long i, Long const& x) {
                                              p_pfsum[i] = x; },
                                          Scan::Type::exclusive, Scan::retSum);
    }
    Long const remote_nnz = m_nnz - local_nnz;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE
        (local_nnz < Long(std::numeric_limits<int>::max()) &&
         remote_nnz < Long(std::numeric_limits<int>::max()),
         "SpMatrix: too many local nonzeros for 32-bit row offsets");

    m_csr_local.resize(nlocalrows, local_nnz);
    container_type<Long> remote_gcols;
    container_type<T> remote_mat;
    container_type<int> remote_row_offset;
#ifdef AMREX_USE_MPI
    remote_gcols.resize(remote_nnz);
    remote_mat.resize(remote_nnz);
    remote_row_offset.resize(nlocalrows+1);
#else
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(remote_nnz == 0,
        "SpMatrix: column index outside the column partition");
#endif
    {
        auto const* pmat = m_csr.mat.data();
        auto* pmat_l = m_csr_local.mat.data();
        auto* pcol_l = m_csr_local.col_index.data();
        auto* pmat_r = remote_mat.data();
        auto* pcol_r = remote_gcols.data();
        ParallelForOMP(m_nnz, [=] AMREX_GPU_DEVICE (Long i)
        {
            auto ps = p_pfsum[i];
            auto local = (pcol[i] >= col_begin && pcol[i] < col_end);
            if (local) {
                pmat_l[ps] = pmat[i];
                pcol_l[ps] = int(pcol[i] - col_begin);
            } else if (pmat_r) {
                pmat_r[i-ps] = pmat[i];
                pcol_r[i-ps] = pcol[i];
            }
        });
        auto const noffset = Long(m_csr.row_offset.size());
        auto const* pro = m_csr.row_offset.data();
        auto* pro_l = m_csr_local.row_offset.data();
        auto* pro_r = remote_row_offset.data();
        auto total_nnz = m_nnz;
        ParallelForOMP(noffset, [=] AMREX_GPU_DEVICE (Long i)
        {
            Long ro_l;
            if (i < noffset-1) {
                ro_l = (pro[i] < total_nnz) ? p_pfsum[pro[i]] : local_nnz;
            } else {
                ro_l = local_nnz;
            }
            pro_l[i] = int(ro_l);
            if (pro_r) { pro_r[i] = int(pro[i] - ro_l); }
        });
        Gpu::streamSynchronize();
    }
    m_csr = csr_type{}; // the unsplit form is no longer needed

#ifdef AMREX_USE_MPI
    // Without remote entries the remote block stays empty (no offsets).
    if (remote_nnz > 0) {
        m_csr_remote.mat = std::move(remote_mat);
        m_csr_remote.col_index.resize(remote_nnz);
        m_csr_remote.row_offset = std::move(remote_row_offset);
        m_csr_remote.nnz = remote_nnz;
        trim_remote_rows();
    }
    update_remote_col_index(remote_gcols, m_csr_remote.col_index, true);
#endif

    m_split = true;
}

#ifdef AMREX_USE_MPI

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::trim_remote_rows ()
{
    BL_PROFILE("SpMatrix::trim_remote_rows");
    // In the remote part, it's expected that some rows don't have
    // nonzeros. So we trim them off, and also save a full copy of the
    // row offsets.
    Long old_size = m_csr_remote.row_offset.size();
    m_ri_ltor.resize(old_size-1);
    m_ri_rtol.resize(old_size-1);
    auto* p_ltor = m_ri_ltor.data();
    auto* p_rtol = m_ri_rtol.data();
    container_type<int> trimmed_row_offset(old_size);
    auto const* p_ro = m_csr_remote.row_offset.data();
    auto* p_tro = trimmed_row_offset.data();
    // This is basically std::unique.
    Long new_size = Scan::PrefixSum<Long>(old_size,
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
                                        if (i+1 < old_size) {
                                            if (p_ro[i+1] > p_ro[i]) {
                                                p_rtol[x] = i;
                                                p_ltor[i] = x;
                                            } else {
                                                p_ltor[i] = -1;
                                            }
                                        }
                                    },
                                    Scan::Type::exclusive, Scan::retSum);
    m_ri_rtol.resize(new_size-1);
    trimmed_row_offset.resize(new_size);
#ifdef AMREX_USE_GPU
    m_ri_rtol.shrink_to_fit();
    trimmed_row_offset.shrink_to_fit();
#endif
    m_remote_row_offset = std::move(trimmed_row_offset);
    std::swap(m_csr_remote.row_offset, m_remote_row_offset);
}

template <typename T, template<typename> class Allocator>
template <typename VL, typename VI>
void SpMatrix<T,Allocator>::update_remote_col_index (VL const& gcols, VI& lcols,
                                                     bool in_device_memory)
{
    int const nprocs = ParallelContext::NProcsSub();

    // This function also needs to update m_remote_cols_*.

    m_remote_cols_v.clear();
    m_remote_cols_vv.clear();
    m_remote_cols_vv.resize(nprocs);
#ifdef AMREX_USE_GPU
    m_remote_cols_dv.clear();
#endif

    Long const n = Long(gcols.size());
    lcols.resize(n);
    if (n == 0) { return; }

    amrex::ignore_unused(in_device_memory);

    Vector<Long> h_gcols;
#ifdef AMREX_USE_GPU
    if (in_device_memory) {
        h_gcols.resize(n);
        Gpu::copyAsync(Gpu::deviceToHost, gcols.begin(), gcols.end(), h_gcols.begin());
        Gpu::streamSynchronize();
    } else
#endif
    {
        h_gcols.assign(gcols.begin(), gcols.end());
    }

    m_remote_cols_v = h_gcols;
    amrex::RemoveDuplicates(m_remote_cols_v);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE
        (Long(m_remote_cols_v.size()) < Long(std::numeric_limits<int>::max()),
         "SpMatrix: too many remote columns for 32-bit indices");

#ifdef AMREX_USE_GPU
    m_remote_cols_dv.resize(m_remote_cols_v.size());
    Gpu::copyAsync(Gpu::hostToDevice,
                   m_remote_cols_v.begin(),
                   m_remote_cols_v.end(),
                   m_remote_cols_dv.data());
#endif

    // Note that amrex::RemoveDuplicates sorts the data.
    auto const& cp = this->m_col_partition.dataVector();
    AMREX_ALWAYS_ASSERT(m_remote_cols_v.front() >= cp.front() &&
                        m_remote_cols_v.back() < cp.back());
    auto it = cp.cbegin();
    for (auto c : m_remote_cols_v) {
        it = std::find_if(it, cp.cend(), [&] (auto x) { return x > c; });
        if (it != cp.cend()) {
            int iproc = int(std::distance(cp.cbegin(),it)) - 1;
            m_remote_cols_vv[iproc].push_back(c);
        } else {
            amrex::Abort("SpMatrix::update_remote_col_index: how did this happen?");
        }
    }

    // Now we convert the remote indices from global to local.
    std::map<Long,int> gtol;
    for (Long i = 0, N = Long(m_remote_cols_v.size()); i < N; ++i) {
        gtol[m_remote_cols_v[i]] = int(i);
    }
    Gpu::PinnedVector<int> h_lcols(n);
    for (Long i = 0; i < n; ++i) {
        h_lcols[i] = gtol[h_gcols[i]];
    }

#ifdef AMREX_USE_GPU
    if (in_device_memory) {
        Gpu::copyAsync(Gpu::hostToDevice, h_lcols.begin(), h_lcols.end(), lcols.begin());
        Gpu::streamSynchronize();
    } else
#endif
    {
        std::copy(h_lcols.begin(), h_lcols.end(), lcols.begin());
    }
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::set_num_neighbors ()
{
    if (m_num_neighbors >= 0) { return; }

    int const nprocs = ParallelContext::NProcsSub();
    auto const mpi_int = ParallelDescriptor::Mpi_typemap<int>::type();
    auto const mpi_comm = ParallelContext::CommunicatorSub();

    amrex::Vector<int> connection(nprocs);
    for (int iproc = 0; iproc < nprocs; ++iproc) {
        connection[iproc] = m_remote_cols_vv[iproc].empty() ? 0 : 1;
    }
    amrex::Vector<int> reduce_scatter_counts(nprocs,1);
    m_num_neighbors = 0;
    BL_MPI_REQUIRE(MPI_Reduce_scatter
                   (connection.data(), &m_num_neighbors, reduce_scatter_counts.data(),
                    mpi_int, MPI_SUM, mpi_comm));
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::prepare_comm_mv (AlgPartition const& col_partition)
{
    BL_PROFILE("SpMatrix::prepare_comm_mv");
    if (m_comm_mv.prepared) { return; }

    // This function needs to be safe when nnz is zero.

    this->split_csr(col_partition);

    int const nprocs = ParallelContext::NProcsSub();
    auto const mpi_tag = ParallelDescriptor::SeqNum();
    auto const mpi_long = ParallelDescriptor::Mpi_typemap<Long>::type();
    auto const mpi_comm = ParallelContext::CommunicatorSub();

    if (m_num_neighbors < 0) { set_num_neighbors(); }

    Vector<MPI_Request> mpi_requests;
    mpi_requests.reserve(nprocs);
    for (int iproc = 0; iproc < nprocs; ++iproc) {
        if ( ! m_remote_cols_vv[iproc].empty()) {
            mpi_requests.push_back(MPI_REQUEST_NULL);
            auto const sz = m_remote_cols_vv[iproc].size();
            if (sz > static_cast<Long>(std::numeric_limits<int>::max())) {
                amrex::Abort("SpMatrix::prepare_comm_mv: remote column payload exceeds MPI int count range.");
            }
            auto const msg_count = static_cast<int>(sz);
            // I need to let other processes know what I need from them.
            BL_MPI_REQUIRE(MPI_Isend(m_remote_cols_vv[iproc].data(),
                                     msg_count,
                                     mpi_long, iproc, mpi_tag, mpi_comm,
                                     &(mpi_requests.back())));
            m_comm_mv.recv_from.push_back(iproc);
            m_comm_mv.recv_counts.push_back(msg_count);
        }
    }

    m_comm_mv.total_counts_recv = Long(m_remote_cols_v.size());

    Vector<Vector<Long>> send_indices(m_num_neighbors);
    m_comm_mv.total_counts_send = 0;
    for (int isend = 0; isend < m_num_neighbors; ++isend) {
        MPI_Status mpi_status;
        BL_MPI_REQUIRE(MPI_Probe(MPI_ANY_SOURCE, mpi_tag, mpi_comm, &mpi_status));
        int receiver = mpi_status.MPI_SOURCE;
        int count;
        BL_MPI_REQUIRE(MPI_Get_count(&mpi_status, mpi_long, &count));
        m_comm_mv.send_to.push_back(receiver);
        m_comm_mv.send_counts.push_back(count);
        send_indices[isend].resize(count);
        BL_MPI_REQUIRE(MPI_Recv(send_indices[isend].data(), count, mpi_long,
                                receiver, mpi_tag, mpi_comm, &mpi_status));
        m_comm_mv.total_counts_send += count;
    }

    m_comm_mv.send_indices.resize(m_comm_mv.total_counts_send);
    Gpu::PinnedVector<Long> send_indices_all;
    send_indices_all.reserve(m_comm_mv.total_counts_send);
    for (auto const& vl : send_indices) {
        for (auto x : vl) {
            send_indices_all.push_back(x);
        }
    }
    Gpu::copyAsync(Gpu::hostToDevice, send_indices_all.begin(), send_indices_all.end(),
                   m_comm_mv.send_indices.begin());
    Gpu::streamSynchronize();

    if (! mpi_requests.empty()) {
        Vector<MPI_Status> mpi_statuses(mpi_requests.size());
        BL_MPI_REQUIRE(MPI_Waitall(int(mpi_requests.size()), mpi_requests.data(),
                                   mpi_statuses.data()));
    }

    m_comm_mv.prepared = true;
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::pack_buffer_mv (AlgVector<T,AllocT> const& v)
{
    auto* pdst = m_comm_mv.send_buffer;
    auto* pidx = m_comm_mv.send_indices.data();
    auto const& vv = v.view();
    auto const nsends = Long(m_comm_mv.send_indices.size());
    ParallelForOMP(nsends, [=] AMREX_GPU_DEVICE (Long i)
    {
        pdst[i] = vv(pidx[i]);
    });
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::unpack_buffer_mv (AlgVector<T,AllocT>& v)
{
    auto const& csr = m_csr_remote;
    if (csr.nnz > 0) {
        T const* AMREX_RESTRICT mat = csr.mat.data();
        auto const* AMREX_RESTRICT col = csr.col_index.data();
        auto const* AMREX_RESTRICT row = csr.row_offset.data();

        auto const* rtol = m_ri_rtol.data();

        auto const* AMREX_RESTRICT px = m_comm_mv.recv_buffer;
        auto      * AMREX_RESTRICT py = v.data();

        auto const nrr = Long(csr.row_offset.size())-1;
        ParallelForOMP(nrr, [=] AMREX_GPU_DEVICE (Long i)
        {
            T r = 0;
            for (Long j = row[i]; j < row[i+1]; ++j) {
                r += mat[j] * px[col[j]];
            }
            py[rtol[i]] += r;
        });
    }
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::unpack_buffer_tr (CommTR const& ctr,
                                              AlgPartition const& col_partition)
{
    m_split = true;
    m_col_partition = col_partition;
    m_col_begin = m_col_partition[ParallelDescriptor::MyProc()  ];
    m_col_end   = m_col_partition[ParallelDescriptor::MyProc()+1];

    m_ri_ltor.resize(numLocalRows(), -1);
    m_remote_cols_vv.resize(ParallelDescriptor::NProcs());

    auto nnz = ctr.total_counts_recv[0];
    if (nnz == 0) { return; }

    m_nnz += nnz;
    auto nb = int(ctr.recv_from.size()); // # of blocked CSRs to be merged
    auto total_local_rows = ctr.total_counts_recv[1];

    // Build compressed row index map
    Vector<Long> ri_map(ctr.recv_buffer_idx_map,
                        ctr.recv_buffer_idx_map + total_local_rows);
    RemoveDuplicates(ri_map);
    Long nrows = ri_map.size(); // # of unique rows.

    // Assembled on the host with global columns, then converted to 32-bit
    // remote indices.
#ifdef AMREX_USE_GPU
    CSR<T,Gpu::PinnedVector,int> csrr;
#else
    auto& csrr = m_csr_remote;
#endif
    Vector<Long> gcols(nnz);
    csrr.mat.resize(nnz);
    csrr.col_index.resize(nnz);
    csrr.row_offset.resize(nrows+1);
    csrr.nnz = nnz;

    // Merge blocks in sender-rank order so that rows stay sorted. Each
    // block covers a contiguous global column range.
    Vector<int> order(nb);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&] (int a, int b) {
        return ctr.recv_from[a] < ctr.recv_from[b]; });

    // Count nnz per compressed row
    Vector<int> row_nnz(nrows, 0);
    for (int i : order) {
        auto nrow_i = ctr.recv_counts[i][1];
        Long const* row_offset = ctr.recv_buffer_row_offset
            +                    ctr.recv_buffer_offset[i][2];
        Long const* idx_map    = ctr.recv_buffer_idx_map
            +                    ctr.recv_buffer_offset[i][3];
        AMREX_ASSERT((row_offset[nrow_i] - row_offset[0]) == ctr.recv_counts[i][0]);

        Long p = 0; // index into ri_map
        for (int lr = 0; lr < nrow_i; ++lr) {
            Long const gr = idx_map[lr];
            while (p < nrows && ri_map[p] < gr) { ++p; }
            AMREX_ASSERT(p < nrows && ri_map[p] == gr);
            // p is now compressed row index
            row_nnz[p] += int(row_offset[lr+1] - row_offset[lr]);
        }
    }

    std::exclusive_scan(row_nnz.begin(), row_nnz.end(), csrr.row_offset.begin(), 0);
    csrr.row_offset.back() = int(csrr.nnz);

    auto rowpos = csrr.row_offset; // make a copy to keep track of offset

    for (int i : order) {
        auto nrow_i = ctr.recv_counts[i][1];
        T    const* mat        = ctr.recv_buffer_mat
            +                    ctr.recv_buffer_offset[i][0];
        Long const* col_index  = ctr.recv_buffer_col_index
            +                    ctr.recv_buffer_offset[i][1];
        Long const* row_offset = ctr.recv_buffer_row_offset
            +                    ctr.recv_buffer_offset[i][2];
        Long const* idx_map    = ctr.recv_buffer_idx_map
            +                    ctr.recv_buffer_offset[i][3];

        Long p = 0; // index into ri_map
        for (int lr = 0; lr < nrow_i; ++lr) {
            Long const gr = idx_map[lr];
            while (p < nrows && ri_map[p] < gr) { ++p; }
            AMREX_ASSERT(p < nrows && ri_map[p] == gr);
            // p is now compressed row index
            auto os_src = row_offset[lr] - row_offset[0];
            auto nvals = row_offset[lr+1] - row_offset[lr];
            auto os_dst = rowpos[p];
            std::memcpy(csrr.mat.data()+os_dst, mat+os_src, sizeof(T)*nvals);
            std::memcpy(gcols.data()+os_dst, col_index+os_src, sizeof(Long)*nvals);
            rowpos[p] += int(nvals);
        }
    }

    m_ri_rtol.resize(nrows);
    Gpu::copyAsync(Gpu::hostToDevice, ri_map.begin(), ri_map.end(), m_ri_rtol.begin());
    {
        auto row_begin = m_row_begin;
        auto* AMREX_RESTRICT ltor = m_ri_ltor.data();
        auto* AMREX_RESTRICT rtol = m_ri_rtol.data();
        ParallelForOMP(nrows, [=] AMREX_GPU_DEVICE (Long i) {
            rtol[i] -= row_begin;
            ltor[rtol[i]] = i;
        });
    }

    update_remote_col_index(gcols, csrr.col_index, false);

#ifdef AMREX_USE_GPU
    amrex::duplicateCSR(Gpu::hostToDevice, m_csr_remote, csrr);
    Gpu::streamSynchronize();
#endif
}

template <typename T, template<typename> class Allocator>
auto SpMatrix<T,Allocator>::fetch_remote_rows_mm (SpMatrix<T,Allocator> const& B)
    -> RemoteRowsMM
{
    // this = A, split by B's row partition. m_comm_mv tells which rows of B
    // other ranks need from us and which we need from them.
    this->prepare_comm_mv(B.partition());

    auto const& cm = m_comm_mv;
    auto const nrecvs = int(cm.recv_from.size());
    auto const nsends = int(cm.send_to.size());
    auto const mpi_long = ParallelDescriptor::Mpi_typemap<Long>::type();
    auto const mpi_t = ParallelDescriptor::Mpi_typemap<T>::type();
    auto const mpi_comm = ParallelContext::CommunicatorSub();

    RemoteRowsMM ext;
    ext.nrows = cm.total_counts_recv; // == m_remote_cols_v.size()
    ext.row_offset.resize(ext.nrows+1);

    auto const b0 = B.m_csr_local.const_view();
    auto const b1 = B.remote_full_const_view();
    Long const b_row_begin = B.m_row_begin;
    Long const nsend_rows = cm.total_counts_send;
    auto const* AMREX_RESTRICT send_idx = cm.send_indices.data();

    // Round 1: number of nonzeros in each requested row.
    Gpu::PinnedVector<Long> h_send(nsend_rows+1);
    {
        Gpu::DeviceVector<Long> d_cnt(nsend_rows);
        auto* pcnt = d_cnt.data();
        ParallelForOMP(nsend_rows, [=] AMREX_GPU_DEVICE (Long i)
        {
            Long const lr = send_idx[i] - b_row_begin;
            pcnt[i] = (b0.row_offset[lr+1] - b0.row_offset[lr])
                +     (b1.row_offset[lr+1] - b1.row_offset[lr]);
        });
        Gpu::copyAsync(Gpu::deviceToHost, d_cnt.begin(), d_cnt.end(), h_send.begin());
        Gpu::streamSynchronize();
    }

    Gpu::PinnedVector<Long> h_recv(ext.nrows+1);
    {
        auto const tag = ParallelDescriptor::SeqNum();
        Vector<MPI_Request> reqs;
        reqs.reserve(nrecvs+nsends);
        Long os = 0;
        for (int i = 0; i < nrecvs; ++i) {
            reqs.push_back(MPI_REQUEST_NULL);
            BL_MPI_REQUIRE(MPI_Irecv(h_recv.data()+os, cm.recv_counts[i], mpi_long,
                                     cm.recv_from[i], tag, mpi_comm, &reqs.back()));
            os += cm.recv_counts[i];
        }
        os = 0;
        for (int i = 0; i < nsends; ++i) {
            reqs.push_back(MPI_REQUEST_NULL);
            BL_MPI_REQUIRE(MPI_Isend(h_send.data()+os, cm.send_counts[i], mpi_long,
                                     cm.send_to[i], tag, mpi_comm, &reqs.back()));
            os += cm.send_counts[i];
        }
        if (! reqs.empty()) {
            Vector<MPI_Status> stats(reqs.size());
            BL_MPI_REQUIRE(MPI_Waitall(int(reqs.size()), reqs.data(), stats.data()));
        }
    }

    // Counts to offsets, in place; per-rank totals.
    auto to_offsets = [] (Gpu::PinnedVector<Long>& v, Long n) {
        Long s = 0;
        for (Long i = 0; i < n; ++i) {
            Long const c = v[i];
            v[i] = s;
            s += c;
        }
        v[n] = s;
    };
    to_offsets(h_recv, ext.nrows);
    to_offsets(h_send, nsend_rows);
    ext.nnz = h_recv[ext.nrows];
    Long const send_nnz = h_send[nsend_rows];

    Vector<Long> recv_nnz(nrecvs), send_nnz_v(nsends);
    for (int i = 0, r = 0; i < nrecvs; ++i) {
        recv_nnz[i] = h_recv[r+cm.recv_counts[i]] - h_recv[r];
        r += cm.recv_counts[i];
        if (recv_nnz[i] >= Long(std::numeric_limits<int>::max())) {
            amrex::Abort("SpMatrix::fetch_remote_rows_mm: message exceeds MPI int count range.");
        }
    }
    for (int i = 0, r = 0; i < nsends; ++i) {
        send_nnz_v[i] = h_send[r+cm.send_counts[i]] - h_send[r];
        r += cm.send_counts[i];
        if (send_nnz_v[i] >= Long(std::numeric_limits<int>::max())) {
            amrex::Abort("SpMatrix::fetch_remote_rows_mm: message exceeds MPI int count range.");
        }
    }

    Gpu::copyAsync(Gpu::hostToDevice, h_recv.begin(), h_recv.end(), ext.row_offset.begin());
    Gpu::DeviceVector<Long> d_send_off(nsend_rows+1);
    Gpu::copyAsync(Gpu::hostToDevice, h_send.begin(), h_send.end(), d_send_off.begin());

    // Round 2: column indices and values.
    auto const tag_c = ParallelDescriptor::SeqNum();
    auto const tag_m = ParallelDescriptor::SeqNum();

    Vector<MPI_Request> rreqs, sreqs;
    if (ext.nnz > 0) {
        ext.col_index = (Long*)The_Pinned_Arena()->alloc(sizeof(Long)*ext.nnz);
        ext.mat = (T*)The_Comms_Arena()->alloc(sizeof(T)*ext.nnz);
        rreqs.reserve(2*nrecvs);
        Long os = 0;
        for (int i = 0; i < nrecvs; ++i) {
            if (recv_nnz[i] > 0) {
                auto n = int(recv_nnz[i]);
                rreqs.push_back(MPI_REQUEST_NULL);
                BL_MPI_REQUIRE(MPI_Irecv(ext.col_index+os, n, mpi_long, cm.recv_from[i],
                                         tag_c, mpi_comm, &rreqs.back()));
                rreqs.push_back(MPI_REQUEST_NULL);
                BL_MPI_REQUIRE(MPI_Irecv(ext.mat+os, n, mpi_t, cm.recv_from[i],
                                         tag_m, mpi_comm, &rreqs.back()));
                os += recv_nnz[i];
            }
        }
    }

    Long* send_col = nullptr;
    T* send_mat = nullptr;
    if (send_nnz > 0) {
        send_col = (Long*)The_Comms_Arena()->alloc(sizeof(Long)*send_nnz);
        send_mat = (T*)The_Comms_Arena()->alloc(sizeof(T)*send_nnz);
        auto const* poff = d_send_off.data();
        Long const b_col_begin = B.m_col_begin;
#ifdef AMREX_USE_GPU
        auto const* b_rcols = B.m_remote_cols_dv.data();
#else
        auto const* b_rcols = B.m_remote_cols_v.data();
#endif
        // Merge the diagonal and off-diagonal parts so that each packed
        // row is sorted by global column index.
        ParallelForOMP(nsend_rows, [=] AMREX_GPU_DEVICE (Long i)
        {
            constexpr Long gmax = std::numeric_limits<Long>::max();
            Long const lr = send_idx[i] - b_row_begin;
            Long p = poff[i];
            Long q0 = b0.row_offset[lr];
            Long q1 = b1.row_offset[lr];
            Long const e0 = b0.row_offset[lr+1];
            Long const e1 = b1.row_offset[lr+1];
            while (q0 < e0 || q1 < e1) {
                Long const g0 = (q0 < e0) ? b0.col_index[q0] + b_col_begin : gmax;
                Long const g1 = (q1 < e1) ? b_rcols[b1.col_index[q1]] : gmax;
                if (g0 < g1) {
                    send_col[p] = g0;
                    send_mat[p] = b0.mat[q0];
                    ++q0;
                } else {
                    send_col[p] = g1;
                    send_mat[p] = b1.mat[q1];
                    ++q1;
                }
                ++p;
            }
        });
        Gpu::streamSynchronize();

        sreqs.reserve(2*nsends);
        Long os = 0;
        for (int i = 0; i < nsends; ++i) {
            if (send_nnz_v[i] > 0) {
                auto n = int(send_nnz_v[i]);
                sreqs.push_back(MPI_REQUEST_NULL);
                BL_MPI_REQUIRE(MPI_Isend(send_col+os, n, mpi_long, cm.send_to[i],
                                         tag_c, mpi_comm, &sreqs.back()));
                sreqs.push_back(MPI_REQUEST_NULL);
                BL_MPI_REQUIRE(MPI_Isend(send_mat+os, n, mpi_t, cm.send_to[i],
                                         tag_m, mpi_comm, &sreqs.back()));
                os += send_nnz_v[i];
            }
        }
    }

    if (! rreqs.empty()) {
        Vector<MPI_Status> stats(rreqs.size());
        BL_MPI_REQUIRE(MPI_Waitall(int(rreqs.size()), rreqs.data(), stats.data()));
    }
    if (! sreqs.empty()) {
        Vector<MPI_Status> stats(sreqs.size());
        BL_MPI_REQUIRE(MPI_Waitall(int(sreqs.size()), sreqs.data(), stats.data()));
    }
    Gpu::streamSynchronize();
    if (send_col) { The_Comms_Arena()->free(send_col); }
    if (send_mat) { The_Comms_Arena()->free(send_mat); }

    return ext;
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::expand_remote_row_offset () const
{
    if (! m_remote_row_offset.empty()) { return; }

    AMREX_ASSERT(m_split);

    auto nrows = numLocalRows();
    m_remote_row_offset.resize(nrows+1);
    auto* pro_full = m_remote_row_offset.data();
    auto const* pro_comp = m_csr_remote.row_offset.data();
    auto const* ri_ltor = m_ri_ltor.data();
    auto nnz_r = m_csr_remote.nnz;
    if (nnz_r == 0 || m_ri_ltor.empty()) {
        // split_csr only builds m_ri_ltor when remote entries exist.
        ParallelForOMP(nrows+1, [=] AMREX_GPU_DEVICE (Long i) { pro_full[i] = 0; });
        Gpu::streamSynchronize();
    } else {
        Scan::PrefixSum<int>(nrows,
                             [=] AMREX_GPU_DEVICE (Long i) -> int {
                                 Long rrow = ri_ltor[i];
                                 if (rrow == -1) {
                                     return int(0);
                                 } else {
                                     return int(pro_comp[rrow+1]-pro_comp[rrow]);
                                 }},
                             [=] AMREX_GPU_DEVICE (Long i, int x) {
                                 if (i == 0) {
                                     pro_full[0] = 0;
                                 }
                                 pro_full[i+1] = x;
                             },
                             Scan::Type::inclusive, Scan::noRetSum);
    }
}

template <typename T, template<typename> class Allocator>
CsrView<T const,int> SpMatrix<T,Allocator>::remote_full_const_view () const
{
    if (m_remote_row_offset.empty()) {
        expand_remote_row_offset();
    }
    auto csr_view = m_csr_remote.const_view();
    csr_view.row_offset = m_remote_row_offset.data();
    csr_view.nrows = numLocalRows();
    return csr_view;
}

#endif

}

#endif
