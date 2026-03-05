#ifndef AMREX_SP_MATRIX_H_
#define AMREX_SP_MATRIX_H_
#include <AMReX_Config.H>

#include <AMReX_AlgPartition.H>
#include <AMReX_AlgVector.H>
#include <AMReX_CSR.H>
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
#include <string>
#include <type_traits>
#include <unordered_map>

namespace amrex {

//! GPU-ready non-owning CSR data container
template <typename T>
struct ParCsr {
    CsrView<T> csr0; // diagonal part
    CsrView<T> csr1; // off-diagonal part
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

template <typename T, template<typename> class Allocator = DefaultAllocator>
class SpMatrix
{
public:
    using value_type = T;
    template <class U> using allocator_type = Allocator<U>;
    template <class U> using container_type = PODVector<U,Allocator<U> >;
    using csr_type = CSR<T,container_type>;
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
     * have duplicagted entries.
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
     * process. It must be valid. Moreover, it should not have duplicagted
     * entries.
     *
     * \param partition Parallel partition information.
     * \param csr CSR represents of the matrix local to this process.
     * \param is_sorted Indicates whether the CSR entries are sorted (within each row).
     */
    void define (AlgPartition partition, csr_type csr, CsrSorted is_sorted);

    [[nodiscard]] AlgPartition const& partition () const { return m_partition; }

    /**
     * \brief Return the column partition used for matrix-vector and
     * matrix-matrix multiplications.
     *
     * This partition corresponds to the row partition of the vector or
     * matrix appearing on the right-hand side of a multiplication.
     */
    [[nodiscard]] AlgPartition const& columnPartition () const { return m_col_partition; }

    [[nodiscard]] Long numLocalRows () const { return m_row_end - m_row_begin; }
    [[nodiscard]] Long numGlobalRows () const { return m_partition.numGlobalRows(); }
    [[nodiscard]] Long numLocalNonZeros () const { return m_nnz; }

    //! Inclusive global index begin.
    [[nodiscard]] Long globalRowBegin () const { return m_row_begin; }
    //! Exclusive global index end.
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
     * for the row specified by the global row index.
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

    //! Return diagonal elements in a square matrix.
    [[nodiscard]] AlgVector<T,AllocT> const& diagonalVector () const;

    //! Return row-sum vector.
    [[nodiscard]] AlgVector<T,AllocT> rowSum () const;

    [[nodiscard]] ParCsr<T      >       parcsr ()      ;
    [[nodiscard]] ParCsr<T const>       parcsr () const;
    [[nodiscard]] ParCsr<T const> const_parcsr () const;

    template <typename U, template<typename> class M, typename N> friend
    void SpMV(AlgVector<U,N>& y, SpMatrix<U,M> const& A, AlgVector<U,N> const& x);

    template <typename U, template<typename> class M> friend
    SpMatrix<U,M> transpose (SpMatrix<U,M> const& A, AlgPartition col_partition);

    template <typename U> friend class AMG;

    //! Private function, but public for cuda
    void define_doit (int nnz_per_row);

    //! Private function, but public for cuda
    template <typename I>
    void define_and_filter_doit (T const* mat, Long const* col_index,
                                 Long nentries, Long const* row_offset);

    void startComm_mv (AlgVector<T,AllocT> const& x);
    void finishComm_mv (AlgVector<T,AllocT>& y);

    void startComm_tr (AlgPartition const& col_partition);
    void finishComm_tr (SpMatrix<T,Allocator>& AT);

private:

    void set_num_neighbors ();

    AlgPartition m_partition;
    AlgPartition m_col_partition;
    Long m_row_begin = 0;
    Long m_row_end = 0;
    Long m_col_begin = 0; //! for diagonal part of mat-vec & mat-mat operations
    Long m_col_end = 0;   //! for diagonal part of mat-vec & mat-mat operations
    Long m_nnz = 0; //! total nnz including both diagonal and off-diagonal
    csr_type m_csr;

    mutable AlgVector<T,AllocT> m_diagonal;

    bool m_split = false; // Has the matrix been split into diagonal and off-diagonal parts?

#ifdef AMREX_USE_MPI
    csr_type m_csr_remote;

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
    container_type<Long> m_ri_ltor; // size: m_csr.nrows()

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

    int m_num_neighbors = -1; // No. of other processes involved in communication

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

    void split_csr (AlgPartition const& col_partition);
    template <typename C>
    void update_remote_col_index (C& csrr, bool in_device_memory);
    void prepare_comm_mv (AlgPartition const& col_partition);
    void pack_buffer_mv (AlgVector<T,AllocT> const& v);
    void unpack_buffer_mv (AlgVector<T,AllocT>& v);
    void unpack_buffer_tr (CommTR const& ctr, AlgPartition const& col_partition);

#endif
};

namespace detail {
template <typename T>
void transpose (CsrView<T> const& csrt, CsrView<T const> const& csr)
{
    Long nrows = csr.nrows;
    Long ncols = csrt.nrows;
    Long nnz = csr.nnz;

    if (nrows <= 0 || ncols <= 0 || nnz <= 0) {
        if (ncols > 0) {
            auto* p = csrt.row_offset;
            ParallelForOMP(ncols+1, [=] AMREX_GPU_DEVICE (Long i) { p[i] = 0; });
        }
        return;
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

    auto* csr_col_index = (int*)The_Arena()->alloc(csr.nnz*sizeof(int));
    auto* csr_row_offset = (int*)The_Arena()->alloc((csr.nrows+1)*sizeof(int));
    auto* csrt_col_index = (int*)The_Arena()->alloc(csrt.nnz*sizeof(int));
    auto* csrt_row_offset = (int*)The_Arena()->alloc((csrt.nrows+1)*sizeof(int));

    ParallelFor(std::max(csr.nnz, csr.nrows+1), [=] AMREX_GPU_DEVICE (Long i)
    {
        if (i < csr.nnz) {
            csr_col_index[i] = int(csr.col_index[i]);
        }
        if (i < csr.nrows+1) {
            csr_row_offset[i] = int(csr.row_offset[i]);
        }
    });

    std::size_t buffer_size;
    AMREX_CUSPARSE_SAFE_CALL(
        cusparseCsr2cscEx2_bufferSize(handle, int(nrows), int(ncols), int(nnz),
                                      csr.mat, csr_row_offset, csr_col_index,
                                      csrt.mat, csrt_row_offset, csrt_col_index,
                                      data_type, CUSPARSE_ACTION_NUMERIC,
                                      CUSPARSE_INDEX_BASE_ZERO,
                                      CUSPARSE_CSR2CSC_ALG1,
                                      &buffer_size));

    auto* pbuffer = (void*)The_Arena()->alloc(buffer_size);

    AMREX_CUSPARSE_SAFE_CALL(
        cusparseCsr2cscEx2(handle, int(nrows), int(ncols), int(nnz),
                           csr.mat, csr_row_offset, csr_col_index,
                           csrt.mat, csrt_row_offset, csrt_col_index,
                           data_type, CUSPARSE_ACTION_NUMERIC,
                           CUSPARSE_INDEX_BASE_ZERO,
                           CUSPARSE_CSR2CSC_ALG1,
                           pbuffer));

    ParallelFor(std::max(csrt.nnz, csrt.nrows+1), [=] AMREX_GPU_DEVICE (Long i)
    {
        if (i < csrt.nnz) {
            csrt.col_index[i] = csrt_col_index[i];
        }
        if (i < csrt.nrows+1) {
            csrt.row_offset[i] = csrt_row_offset[i];
        }
    });

    Gpu::streamSynchronize();
    AMREX_CUSPARSE_SAFE_CALL(cusparseDestroy(handle));
    The_Arena()->free(pbuffer);
    The_Arena()->free(csr_row_offset);
    The_Arena()->free(csr_col_index);
    The_Arena()->free(csrt_row_offset);
    The_Arena()->free(csrt_col_index);

#elif defined(AMREX_USE_HIP)

    rocsparse_handle handle;
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_create_handle(&handle));
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_set_stream(handle, Gpu::gpuStream()));

    AMREX_ALWAYS_ASSERT(nrows < Long(std::numeric_limits<rocsparse_int>::max()) &&
                        ncols < Long(std::numeric_limits<rocsparse_int>::max()) &&
                        nnz   < Long(std::numeric_limits<rocsparse_int>::max()));

    rocsparse_int* csr_col_index;
    rocsparse_int* csr_row_offset;
    rocsparse_int* csrt_col_index;
    rocsparse_int* csrt_row_offset;
    if (std::is_same_v<rocsparse_int,Long>) {
        csr_col_index = (rocsparse_int*)csr.col_index;
        csr_row_offset = (rocsparse_int*)csr.row_offset;
        csrt_col_index = (rocsparse_int*)csrt.col_index;
        csrt_row_offset = (rocsparse_int*)csrt.row_offset;
    } else {
        csr_col_index = (rocsparse_int*)The_Arena()->alloc(csr.nnz*sizeof(rocsparse_int));
        csr_row_offset = (rocsparse_int*)The_Arena()->alloc((csr.nrows+1)*sizeof(rocsparse_int));
        csrt_col_index = (rocsparse_int*)The_Arena()->alloc(csrt.nnz*sizeof(rocsparse_int));
        csrt_row_offset = (rocsparse_int*)The_Arena()->alloc((csrt.nrows+1)*sizeof(rocsparse_int));
        ParallelFor(std::max(csr.nnz, csr.nrows+1), [=] AMREX_GPU_DEVICE (Long i)
        {
            if (i < csr.nnz) {
                csr_col_index[i] = rocsparse_int(csr.col_index[i]);
            }
            if (i < csr.nrows+1) {
                csr_row_offset[i] = rocsparse_int(csr.row_offset[i]);
            }
        });
    }

    std::size_t buffer_size;
    AMREX_ROCSPARSE_SAFE_CALL(
        rocsparse_csr2csc_buffer_size(handle, rocsparse_int(nrows),
                                      rocsparse_int(ncols), rocsparse_int(nnz),
                                      csr_row_offset, csr_col_index,
                                      rocsparse_action_numeric,
                                      &buffer_size));

    auto* pbuffer = (void*)The_Arena()->alloc(buffer_size);

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

    if (! std::is_same_v<rocsparse_int,Long>) {
        ParallelFor(std::max(csrt.nnz, csrt.nrows+1), [=] AMREX_GPU_DEVICE (Long i)
        {
            if (i < csrt.nnz) {
                csrt.col_index[i] = csrt_col_index[i];
            }
            if (i < csrt.nrows+1) {
                csrt.row_offset[i] = csrt_row_offset[i];
            }
        });
    }

    Gpu::streamSynchronize();
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_destroy_handle(handle));
    The_Arena()->free(pbuffer);
    if (! std::is_same_v<rocsparse_int,Long>) {
        The_Arena()->free(csr_row_offset);
        The_Arena()->free(csr_col_index);
        The_Arena()->free(csrt_row_offset);
        The_Arena()->free(csrt_col_index);
    }

#elif defined(AMREX_USE_SYCL)

    mkl::sparse::matrix_handle_t handle_in{};
    mkl::sparse::matrix_handle_t handle_out{};
    mkl::sparse::init_matrix_handle(&handle_in);
    mkl::sparse::init_matrix_handle(&handle_out);

#if defined(INTEL_MKL_VERSION) && (INTEL_MKL_VERSION < 20250300)
    amrex::ignore_unused(nnz);
    mkl::sparse::set_csr_data(Gpu::Device::streamQueue(), handle_in, nrows, ncols,
                              mkl::index_base::zero, (Long*)csr.row_offset,
                              (Long*)csr.col_index, (T*)csr.mat);
    mkl::sparse::set_csr_data(Gpu::Device::streamQueue(), handle_out, ncols, nrows,
                              mkl::index_base::zero, (Long*)csrt.row_offset,
                              (Long*)csrt.col_index, (T*)csrt.mat);
#else
    mkl::sparse::set_csr_data(Gpu::Device::streamQueue(), handle_in, nrows, ncols, nnz,
                              mkl::index_base::zero, (Long*)csr.row_offset,
                              (Long*)csr.col_index, (T*)csr.mat);
    mkl::sparse::set_csr_data(Gpu::Device::streamQueue(), handle_out, ncols, nrows, nnz,
                              mkl::index_base::zero, (Long*)csrt.row_offset,
                              (Long*)csrt.col_index, (T*)csrt.mat);
#endif

    mkl::sparse::omatcopy(Gpu::Device::streamQueue(), mkl::transpose::trans,
                          handle_in, handle_out);

    mkl::sparse::release_matrix_handle(Gpu::Device::streamQueue(), &handle_in);
    auto ev = mkl::sparse::release_matrix_handle(Gpu::Device::streamQueue(), &handle_out);
    ev.wait();

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
            csrt.col_index[dest] = i;
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
    define_doit(nnz_per_row);
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::define (AlgPartition partition, csr_type csr,
                                    CsrSorted is_sorted)
{
    m_partition = std::move(partition);
    m_row_begin = m_partition[ParallelDescriptor::MyProc()];
    m_row_end = m_partition[ParallelDescriptor::MyProc()+1];
    m_nnz = csr.nnz;
    m_csr = std::move(csr);
    if (! is_sorted) { m_csr.sort(); }
}

template <typename T, template<typename> class Allocator>
void
SpMatrix<T,Allocator>::define_doit (int nnz_per_row)
{
    if (nnz_per_row <= 0) { return; };

    AMREX_ALWAYS_ASSERT(m_split == false);

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
    m_csr.sort();
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
    ParallelFor(std::max(nentries,nlocalrows), [=] AMREX_GPU_DEVICE (Long i)
    {
        if (i < nentries) {
            if (col_index[i] >= 0 && mat[i] != 0) {
                pmat[ps[i]] = mat[i];
                pcol[ps[i]] = col_index[i];
            }
        }
        if (i < nlocalrows) {
            prow[i] = ps[row_offset[i]];
            if (i == nlocalrows - 1) {
                prow[nlocalrows] = actual_nnz;
            }
        }
    });
    Gpu::streamSynchronize();
}

template <typename T, template<typename> class Allocator>
void
SpMatrix<T,Allocator>::printToFile (std::string const& file) const
{
#ifdef AMREX_USE_GPU
    CSR<T,Gpu::PinnedVector> csr;
    amrex::duplicateCSR(Gpu::deviceToHost, csr, m_csr);

#  ifdef AMREX_USE_MPI
    CSR<T,Gpu::PinnedVector> csr_r;
    if (m_split) {
        amrex::duplicateCSR(Gpu::deviceToHost, csr_r, m_csr_remote);
    }

    Gpu::PinnedVector<Long> ri_ltor(m_ri_ltor.size());
    if (m_split) {
        Gpu::copyAsync(Gpu::deviceToHost,
                       m_ri_ltor.begin(),
                       m_ri_ltor.end(),
                       ri_ltor.begin());
    }
    auto const& remote_cols = m_remote_cols_v;
#  endif

    Gpu::streamSynchronize();

#else

    auto const& csr = m_csr;
#  ifdef AMREX_USE_MPI
    auto const& csr_r = m_csr_remote;
    auto const& ri_ltor = m_ri_ltor;
    auto const& remote_cols = m_remote_cols_v;
#  endif

#endif

    Long nnz = m_csr.nnz;
#ifdef AMREX_USE_MPI
    nnz += m_csr_remote.nnz;
#endif

    std::ofstream ofs(file+"."+std::to_string(ParallelDescriptor::MyProc()));
    ofs << m_row_begin << " " << m_row_end << " " << nnz << "\n";
    for (Long i = 0, nrows = numLocalRows(); i < nrows; ++i) {
        Long nnz_row    = csr.row_offset[i+1]  - csr.row_offset[i];
        T    const* mat = csr.mat.data()       + csr.row_offset[i];
        Long const* col = csr.col_index.data() + csr.row_offset[i];
        for (Long j = 0; j < nnz_row; ++j) {
            ofs << i+m_row_begin << " " << col[j]+m_col_begin << " " << mat[j] << "\n";
        }
#ifdef AMREX_USE_MPI
        if (i < Long(ri_ltor.size()) && ri_ltor[i] >= 0) {
            Long ii = ri_ltor[i];
            nnz_row = csr_r.row_offset[ii+1]  - csr_r.row_offset[ii];
            mat     = csr_r.mat.data()        + csr_r.row_offset[ii];
            col     = csr_r.col_index.data()  + csr_r.row_offset[ii];
            for (Long j = 0; j < nnz_row; ++j) {
                ofs << i+m_row_begin << " " << remote_cols[col[j]] << " " << mat[j] << "\n";
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
        m_diagonal.define(this->partition());
        auto* AMREX_RESTRICT p = m_diagonal.data();
        auto const* AMREX_RESTRICT mat = m_csr.mat.data();
        auto const* AMREX_RESTRICT col = m_csr.col_index.data();
        auto const* AMREX_RESTRICT row = m_csr.row_offset.data();
        auto offset = m_split ? Long(0) : m_row_begin; // assuming square matrix
        Long nrows = this->numLocalRows();
        ParallelForOMP(nrows, [=] AMREX_GPU_DEVICE (Long i)
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
AlgVector<T,Allocator<T>> SpMatrix<T,Allocator>::rowSum () const
{
    AlgVector<T,Allocator<T>> r(this->partition());
    auto* p = r.data();
    auto const& a = this->const_parcsr();
    ParallelForOMP(r.numLocalRows(), [=] AMREX_GPU_DEVICE (Long i)
    {
        T s = 0;
        for (auto idx = a.csr0.row_offset[i];
                  idx < a.csr0.row_offset[i+1]; ++idx) {
            s += a.csr0.mat[idx];
        }
        if (a.csr1.nnz > 0) {
            for (auto idx = a.csr1.row_offset[i];
                      idx < a.csr1.row_offset[i+1]; ++idx) {
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
    return ParCsr<T>{m_csr.view(),
#ifdef AMREX_USE_MPI
                     m_csr_remote.view(),
#else
                     CsrView<T>{},
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
    return ParCsr<U>{m_csr.const_view(),
#ifdef AMREX_USE_MPI
                     m_csr_remote.const_view(),
#else
                     CsrView<U>{},
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
#ifndef AMREX_USE_MPI
    amrex::ignore_unused(x);
#else
    if (this->partition().numActiveProcs() <= 1) { return; }

    this->prepare_comm_mv(x.partition());

    auto const mpi_tag = ParallelDescriptor::SeqNum();
    auto const mpi_t_type = ParallelDescriptor::Mpi_typemap<T>::type(); // NOLINT(readability-qualified-auto)
    auto const mpi_comm = ParallelContext::CommunicatorSub(); // NOLINT(readability-qualified-auto)

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
    if (this->numLocalRows() == 0) { return; }

#ifndef AMREX_USE_MPI
    amrex::ignore_unused(y);
#else
    if (this->partition().numActiveProcs() <= 1) { return; }

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
void SpMatrix<T,Allocator>::startComm_tr (AlgPartition const& col_partition)
{
#ifdef AMREX_USE_MPI
    if (this->partition().numActiveProcs() <= 1) { return; }

    this->split_csr(col_partition);

    int const nprocs = ParallelContext::NProcsSub();
    auto const mpi_tag = ParallelDescriptor::SeqNum();
    auto const mpi_long = ParallelDescriptor::Mpi_typemap<Long>::type(); // NOLINT(readability-qualified-auto)
    auto const mpi_t = ParallelDescriptor::Mpi_typemap<T>::type(); // NOLINT(readability-qualified-auto)
    auto const mpi_comm = ParallelContext::CommunicatorSub(); // NOLINT(readability-qualified-auto)

    // transpose the off-diagonal part
    if (m_csr_remote.nnz > 0) {
        m_comm_tr.csrt.nnz = m_csr_remote.nnz;
        m_comm_tr.csrt.nrows = m_remote_cols_v.size();
        m_comm_tr.csrt.mat = (T*)The_Comms_Arena()->alloc
            (sizeof(T)*m_comm_tr.csrt.nnz);
        m_comm_tr.csrt.col_index = (Long*)The_Comms_Arena()->alloc
            (sizeof(Long)*m_comm_tr.csrt.nnz);
        m_comm_tr.csrt.row_offset = (Long*)The_Comms_Arena()->alloc
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
#ifdef AMREX_USE_MPI
    if (this->partition().numActiveProcs() <= 1) { return; }

    if (! m_comm_tr.recv_reqs.empty()) {
        Vector<MPI_Status> mpi_statuses(m_comm_tr.recv_reqs.size());
        BL_MPI_REQUIRE(MPI_Waitall(int(m_comm_tr.recv_reqs.size()),
                                   m_comm_tr.recv_reqs.data(),
                                   mpi_statuses.data()));
    }

    AT.unpack_buffer_tr(m_comm_tr, this->m_partition);

    if (! m_comm_tr.send_reqs.empty()) {
        Vector<MPI_Status> mpi_statuses(m_comm_tr.send_reqs.size());
        BL_MPI_REQUIRE(MPI_Waitall(int(m_comm_tr.send_reqs.size()),
                                   m_comm_tr.send_reqs.data(),
                                   mpi_statuses.data()));
    }

    if (m_comm_tr.csrt.nnz > 0) {
        The_Comms_Arena()->free(m_comm_tr.csrt.mat);
        The_Comms_Arena()->free(m_comm_tr.csrt.col_index);
        The_Comms_Arena()->free(m_comm_tr.csrt.row_offset);
    }
    if (m_comm_tr.recv_buffer_mat) {
        The_Pinned_Arena()->free(m_comm_tr.recv_buffer_mat);
        The_Pinned_Arena()->free(m_comm_tr.recv_buffer_col_index);
        The_Pinned_Arena()->free(m_comm_tr.recv_buffer_row_offset);
        The_Pinned_Arena()->free(m_comm_tr.recv_buffer_idx_map);
    }
    m_comm_tr = CommTR{};
#else
    amrex::ignore_unused(AT);
#endif
}

#ifdef AMREX_USE_MPI

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::split_csr (AlgPartition const& col_partition)
{
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

    // This function needs to be safe when nnz is zero.

    // We need to split the matrix into two parts, a diagonal part for pure
    // local operations and another part for remote operations in
    // matrix-vector or matrix-matrix multiplication.

    Long local_nnz;
    Gpu::DeviceVector<Long> pfsum(m_nnz);
    auto* p_pfsum = pfsum.data();
    auto col_begin = m_col_begin;
    auto col_end = m_col_end;
    if (m_csr.nnz < Long(std::numeric_limits<int>::max())) {
        auto const* pcol = m_csr.col_index.data();
        local_nnz = Scan::PrefixSum<int>(int(m_nnz),
                                         [=] AMREX_GPU_DEVICE (int i) -> int {
                                             return (pcol[i] >= col_begin &&
                                                     pcol[i] <  col_end); },
                                         [=] AMREX_GPU_DEVICE (int i, int const& x) {
                                             p_pfsum[i] = x; },
                                         Scan::Type::exclusive, Scan::retSum);
    } else {
        auto const* pcol = m_csr.col_index.data();
        local_nnz = Scan::PrefixSum<Long>(m_nnz,
                                          [=] AMREX_GPU_DEVICE (Long i) -> Long {
                                              return (pcol[i] >= col_begin &&
                                                      pcol[i] <  col_end); },
                                          [=] AMREX_GPU_DEVICE (Long i, Long const& x) {
                                              p_pfsum[i] = x; },
                                          Scan::Type::exclusive, Scan::retSum);
    }

    m_csr.nnz = local_nnz;
    Long remote_nnz = m_nnz - local_nnz;
    m_csr_remote.nnz = remote_nnz;

    if (local_nnz != m_nnz) {
        m_csr_remote.mat.resize(remote_nnz);
        m_csr_remote.col_index.resize(remote_nnz);
        container_type<T> new_mat(local_nnz);
        container_type<Long> new_col(local_nnz);
        auto const* pmat = m_csr.mat.data();
        auto const* pcol = m_csr.col_index.data();
        auto* pmat_l = new_mat.data();
        auto* pcol_l = new_col.data();
        auto* pmat_r = m_csr_remote.mat.data();
        auto* pcol_r = m_csr_remote.col_index.data();
        ParallelForOMP(m_nnz, [=] AMREX_GPU_DEVICE (Long i)
        {
            auto ps = p_pfsum[i];
            auto local = (pcol[i] >= col_begin &&
                          pcol[i] <  col_end);
            if (local) {
                pmat_l[ps] = pmat[i];
                pcol_l[ps] = pcol[i] - col_begin; // shift the column index to local
            } else {
                pmat_r[i-ps] = pmat[i];
                pcol_r[i-ps] = pcol[i];
            }
        });
        auto noffset = Long(m_csr.row_offset.size());
        auto* pro = m_csr.row_offset.data();
        m_csr_remote.row_offset.resize(noffset);
        auto* pro_r = m_csr_remote.row_offset.data();
        ParallelForOMP(noffset, [=] AMREX_GPU_DEVICE (Long i)
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
        m_csr.mat.swap(new_mat);
        m_csr.col_index.swap(new_col);

        // In the remote part, it's expected that some rows don't have
        // nonzeros. So we trim them off.
        {
            Long old_size = m_csr_remote.row_offset.size();
            m_ri_ltor.resize(old_size-1);
            m_ri_rtol.resize(old_size-1);
            auto* p_ltor = m_ri_ltor.data();
            auto* p_rtol = m_ri_rtol.data();
            container_type<Long> trimmed_row_offset(old_size);
            auto const* p_ro = m_csr_remote.row_offset.data();
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
            }

            m_ri_rtol.resize(new_size-1);
            trimmed_row_offset.resize(new_size);
#ifdef AMREX_USE_GPU
            m_ri_rtol.shrink_to_fit();
            trimmed_row_offset.shrink_to_fit();
#endif
            m_csr_remote.row_offset.swap(trimmed_row_offset);
        }

    } else if (col_begin > 0) {
        auto* pcol = m_csr.col_index.data();
        ParallelForOMP(m_nnz, [=] AMREX_GPU_DEVICE (Long i) { pcol[i] -= col_begin; });
    }

    update_remote_col_index(m_csr_remote, true);

    m_split = true;
}

template <typename T, template<typename> class Allocator>
template <typename C>
void SpMatrix<T,Allocator>::update_remote_col_index (C& csrr, bool in_device_memory)
{
    int const nprocs = ParallelContext::NProcsSub();

    // This function also needs to update m_remote_cols_*.

    m_remote_cols_v.clear();
    m_remote_cols_vv.clear();
    m_remote_cols_vv.resize(nprocs);
#ifdef AMREX_USE_GPU
    m_remote_cols_dv.clear();
#endif

    if (csrr.nnz == 0) { return; }

    amrex::ignore_unused(in_device_memory);

#ifdef AMREX_USE_GPU
    if (in_device_memory) {
        m_remote_cols_v.resize(csrr.nnz);
        Gpu::copyAsync(Gpu::deviceToHost,
                       csrr.col_index.begin(),
                       csrr.col_index.end(),
                       m_remote_cols_v.begin());
        Gpu::streamSynchronize();
    } else
#endif
    {
        m_remote_cols_v.assign(csrr.col_index.begin(),
                               csrr.col_index.end());
    }

    amrex::RemoveDuplicates(m_remote_cols_v);

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
    std::map<Long,Long> gtol;
    for (Long i = 0, N = Long(m_remote_cols_v.size()); i < N; ++i) {
        gtol[m_remote_cols_v[i]] = i;
    }

#ifdef AMREX_USE_GPU
    if (in_device_memory) {
        Gpu::PinnedVector<Long> host_col_index(csrr.nnz);
        Gpu::copyAsync(Gpu::deviceToHost,
                       csrr.col_index.begin(),
                       csrr.col_index.end(),
                       host_col_index.begin());
        Gpu::streamSynchronize();
        for (auto& c : host_col_index) {
            c = gtol[c];
        }
        Gpu::copyAsync(Gpu::hostToDevice,
                       host_col_index.begin(),
                       host_col_index.end(),
                       csrr.col_index.begin());
        Gpu::streamSynchronize();
    } else
#endif
    {
        for (auto& c : csrr.col_index) {
            c = gtol[c];
        }
    }
}

template <typename T, template<typename> class Allocator>
void SpMatrix<T,Allocator>::set_num_neighbors ()
{
    if (m_num_neighbors >= 0) { return; }

    int const nprocs = ParallelContext::NProcsSub();
    auto const mpi_int = ParallelDescriptor::Mpi_typemap<int>::type(); // NOLINT(readability-qualified-auto)
    auto const mpi_comm = ParallelContext::CommunicatorSub(); // NOLINT(readability-qualified-auto)

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
    if (m_comm_mv.prepared) { return; }

    // This function needs to be safe when nnz is zero.

    this->split_csr(col_partition);

    int const nprocs = ParallelContext::NProcsSub();
    auto const mpi_tag = ParallelDescriptor::SeqNum();
    auto const mpi_long = ParallelDescriptor::Mpi_typemap<Long>::type(); // NOLINT(readability-qualified-auto)
    auto const mpi_comm = ParallelContext::CommunicatorSub(); // NOLINT(readability-qualified-auto)

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

    m_ri_ltor.resize(m_csr.nrows(), -1);
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

#ifdef AMREX_USE_GPU
    CSR<T,Gpu::PinnedVector> csrr;
#else
    auto& csrr = m_csr_remote;
#endif
    csrr.mat.resize(nnz);
    csrr.col_index.resize(nnz);
    csrr.row_offset.resize(nrows+1);
    csrr.nnz = nnz;

    // Count nnz per compressed row
    Vector<int> row_nnz(nrows, 0);
    for (int i = 0; i < nb; ++i) {
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

    csrr.row_offset[0] = 0;
    std::partial_sum(row_nnz.begin(), row_nnz.end(), csrr.row_offset.begin()+1);
    AMREX_ASSERT(csrr.nnz == csrr.row_offset.back());

    auto rowpos = csrr.row_offset; // make a copy to keep track of offset

    for (int i = 0; i < nb; ++i) {
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
            std::memcpy(csrr.      mat.data()+os_dst,       mat+os_src,
                        sizeof(T)   *nvals);
            std::memcpy(csrr.col_index.data()+os_dst, col_index+os_src,
                        sizeof(Long)*nvals);

            rowpos[p] += nvals;
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

    // The column index in csrr is still global.
    update_remote_col_index(csrr, false);

#ifdef AMREX_USE_GPU
    amrex::duplicateCSR(Gpu::hostToDevice, m_csr_remote, csrr);
    Gpu::streamSynchronize();
#endif
}

#endif

}

#endif
