#ifndef AMREX_SPGEMM_H_
#define AMREX_SPGEMM_H_
#include <AMReX_Config.H>

#include <AMReX_Algorithm.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Exception.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_OpenMP.H>
#include <AMReX_SpMatrix.H>

#include <algorithm>
#include <limits>
#include <utility>

namespace amrex::detail {

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

// Gustavson's algorithm. Pass 1 counts, pass 2 accumulates into a dense
// per-thread accumulator. Output rows are sorted.
template <typename T, template<typename> class V, typename I>
CSR<T,V,I> spgemm_local_cpu (Long nrows, Long ncols,
                             CsrView<T const,I> const& A, CsrView<T const,I> const& B)
{
    CSR<T,V,I> C;
    C.row_offset.resize(nrows+1);
    I* AMREX_RESTRICT crow = C.row_offset.data();
    AMREX_ASSUME(crow != nullptr); // gcc -Wnull-dereference false positive
    crow[0] = 0;

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    {
        Vector<Long> marker(ncols, -1);
#ifdef AMREX_USE_OMP
#pragma omp for schedule(dynamic,64)
#endif
        for (Long i = 0; i < nrows; ++i) {
            Long cnt = 0;
            for (Long ap = A.row_offset[i]; ap < A.row_offset[i+1]; ++ap) {
                Long const k = A.col_index[ap];
                for (Long bp = B.row_offset[k]; bp < B.row_offset[k+1]; ++bp) {
                    Long const j = B.col_index[bp];
                    if (marker[j] != i) {
                        marker[j] = i;
                        ++cnt;
                    }
                }
            }
            crow[i+1] = I(cnt);
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
    I*    AMREX_RESTRICT ccol = C.col_index.data();
    T*    AMREX_RESTRICT cmat = C.mat.data();

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    {
        Vector<Long> marker(ncols, -1);
        Vector<T> acc(ncols);
        Vector<I> cols;
#ifdef AMREX_USE_OMP
#pragma omp for schedule(dynamic,64)
#endif
        for (Long i = 0; i < nrows; ++i) {
            cols.clear();
            for (Long ap = A.row_offset[i]; ap < A.row_offset[i+1]; ++ap) {
                Long const k = A.col_index[ap];
                T const a = A.mat[ap];
                for (Long bp = B.row_offset[k]; bp < B.row_offset[k+1]; ++bp) {
                    Long const j = B.col_index[bp];
                    if (marker[j] != i) {
                        marker[j] = i;
                        cols.push_back(j);
                        acc[j] = a * B.mat[bp];
                    } else {
                        acc[j] += a * B.mat[bp];
                    }
                }
            }
            std::sort(cols.begin(), cols.end());
            Long p = crow[i];
            for (I j : cols) {
                ccol[p] = j;
                cmat[p] = acc[j];
                ++p;
            }
        }
    }

    return C;
}

#elif defined(AMREX_USE_CUDA)

inline void spgemm_cusparse_check (cusparseStatus_t status)
{
    if (status == CUSPARSE_STATUS_ALLOC_FAILED ||
        status == CUSPARSE_STATUS_INSUFFICIENT_RESOURCES) {
        (void)cudaGetLastError(); // clear a failed internal cudaMalloc
        throw OutOfMemoryError("SpGEMM: cuSPARSE ran out of memory");
    }
    AMREX_CUSPARSE_SAFE_CALL(status);
}

template <typename T, template<typename> class V, typename I>
CSR<T,V,I> spgemm_local_cusparse (Long nrows, Long ncols,
                                  CsrView<T const,I> const& A, CsrView<T const,I> const& B)
{
    static_assert(std::is_same_v<I,int>, "spgemm_local_cusparse: 32-bit indices only");

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

    AMREX_ALWAYS_ASSERT(nrows < Long(std::numeric_limits<int>::max()) &&
                        ncols < Long(std::numeric_limits<int>::max()));

    CSR<T,V,I> C;
    C.row_offset.resize(nrows+1);

    // Releases everything, also when an OutOfMemoryError is thrown.
    struct Resources {
        cusparseHandle_t handle = nullptr;
        cusparseSpMatDescr_t mat_A = nullptr, mat_B = nullptr, mat_C = nullptr;
        cusparseSpGEMMDescr_t descr = nullptr;
        void* buffer1 = nullptr;
        void* buffer2 = nullptr;
        Resources () = default;
        Resources (Resources const&) = delete;
        Resources& operator= (Resources const&) = delete;
        ~Resources () {
            Gpu::streamSynchronize();
            if (descr) { cusparseSpGEMM_destroyDescr(descr); }
            if (mat_A) { cusparseDestroySpMat(mat_A); }
            if (mat_B) { cusparseDestroySpMat(mat_B); }
            if (mat_C) { cusparseDestroySpMat(mat_C); }
            if (handle) { cusparseDestroy(handle); }
            if (buffer1) { The_Arena()->free(buffer1); }
            if (buffer2) { The_Arena()->free(buffer2); }
        }
    } r;

    AMREX_CUSPARSE_SAFE_CALL(cusparseCreate(&r.handle));
    AMREX_CUSPARSE_SAFE_CALL(cusparseSetStream(r.handle, Gpu::gpuStream()));

    constexpr cusparseIndexType_t index_type = CUSPARSE_INDEX_32I;
    void* rowA = (void*)A.row_offset;
    void* colA = (void*)A.col_index;
    void* rowB = (void*)B.row_offset;
    void* colB = (void*)B.col_index;
    void* rowC = (void*)C.row_offset.data();

    AMREX_CUSPARSE_SAFE_CALL
        (cusparseCreateCsr(&r.mat_A, nrows, B.nrows, A.nnz, rowA, colA, (void*)A.mat,
                           index_type, index_type, CUSPARSE_INDEX_BASE_ZERO, data_type));
    AMREX_CUSPARSE_SAFE_CALL
        (cusparseCreateCsr(&r.mat_B, B.nrows, ncols, B.nnz, rowB, colB, (void*)B.mat,
                           index_type, index_type, CUSPARSE_INDEX_BASE_ZERO, data_type));
    AMREX_CUSPARSE_SAFE_CALL
        (cusparseCreateCsr(&r.mat_C, nrows, ncols, 0, rowC, nullptr, nullptr,
                           index_type, index_type, CUSPARSE_INDEX_BASE_ZERO, data_type));

    AMREX_CUSPARSE_SAFE_CALL(cusparseSpGEMM_createDescr(&r.descr));

    T alpha = T(1);
    T beta  = T(0);
    cusparseOperation_t op = CUSPARSE_OPERATION_NON_TRANSPOSE;
    auto const alg = CUSPARSE_SPGEMM_DEFAULT;

    std::size_t buffer_size1 = 0;
    spgemm_cusparse_check
        (cusparseSpGEMM_workEstimation(r.handle, op, op, &alpha, r.mat_A, r.mat_B, &beta,
                                       r.mat_C, data_type, alg, r.descr,
                                       &buffer_size1, nullptr));
    r.buffer1 = (void*)The_Arena()->alloc(buffer_size1);
    spgemm_cusparse_check
        (cusparseSpGEMM_workEstimation(r.handle, op, op, &alpha, r.mat_A, r.mat_B, &beta,
                                       r.mat_C, data_type, alg, r.descr,
                                       &buffer_size1, r.buffer1));

    std::size_t buffer_size2 = 0;
    spgemm_cusparse_check
        (cusparseSpGEMM_compute(r.handle, op, op, &alpha, r.mat_A, r.mat_B, &beta,
                                r.mat_C, data_type, alg, r.descr,
                                &buffer_size2, nullptr));
    r.buffer2 = (void*)The_Arena()->alloc(buffer_size2);
    spgemm_cusparse_check
        (cusparseSpGEMM_compute(r.handle, op, op, &alpha, r.mat_A, r.mat_B, &beta,
                                r.mat_C, data_type, alg, r.descr,
                                &buffer_size2, r.buffer2));

    std::int64_t c_nrows, c_ncols, c_nnz;
    AMREX_CUSPARSE_SAFE_CALL(cusparseSpMatGetSize(r.mat_C, &c_nrows, &c_ncols, &c_nnz));
    AMREX_ALWAYS_ASSERT(c_nrows == nrows && c_ncols == ncols);
    AMREX_ALWAYS_ASSERT(c_nnz < std::int64_t(std::numeric_limits<int>::max()));

    C.mat.resize(c_nnz);
    C.col_index.resize(c_nnz);
    C.nnz = c_nnz;
    void* colC = (void*)C.col_index.data();
    AMREX_CUSPARSE_SAFE_CALL(cusparseCsrSetPointers(r.mat_C, rowC, colC, (void*)C.mat.data()));

    // cuSPARSE guarantees sorted column indices in C.
    AMREX_CUSPARSE_SAFE_CALL
        (cusparseSpGEMM_copy(r.handle, op, op, &alpha, r.mat_A, r.mat_B, &beta, r.mat_C,
                             data_type, alg, r.descr));

    return C;
}

#elif defined(AMREX_USE_HIP)

inline void spgemm_rocsparse_check (rocsparse_status status)
{
    if (status == rocsparse_status_memory_error) {
        (void)hipGetLastError(); // clear a failed internal hipMalloc
        throw OutOfMemoryError("SpGEMM: rocSPARSE ran out of memory");
    }
    AMREX_ROCSPARSE_SAFE_CALL(status);
}

template <typename T, template<typename> class V, typename I>
CSR<T,V,I> spgemm_local_rocsparse (Long nrows, Long ncols,
                                   CsrView<T const,I> const& A, CsrView<T const,I> const& B)
{
    static_assert(std::is_same_v<I,int>, "spgemm_local_rocsparse: 32-bit indices only");

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

    // Releases everything, also when an OutOfMemoryError is thrown.
    struct Resources {
        rocsparse_handle handle = nullptr;
        rocsparse_spmat_descr mat_A = nullptr, mat_B = nullptr, mat_C = nullptr,
            mat_D = nullptr;
        void* buffer = nullptr;
        Resources () = default;
        Resources (Resources const&) = delete;
        Resources& operator= (Resources const&) = delete;
        ~Resources () {
            Gpu::streamSynchronize();
            for (auto m : {mat_A, mat_B, mat_C, mat_D}) {
                if (m) { rocsparse_destroy_spmat_descr(m); }
            }
            if (handle) { rocsparse_destroy_handle(handle); }
            if (buffer) { The_Arena()->free(buffer); }
        }
    } r;

    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_create_handle(&r.handle));
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_set_stream(r.handle, Gpu::gpuStream()));

    AMREX_ROCSPARSE_SAFE_CALL
        (rocsparse_create_csr_descr(&r.mat_A, nrows, B.nrows, A.nnz,
                                    (void*)A.row_offset, (void*)A.col_index, (void*)A.mat,
                                    index_type, index_type, index_base, data_type));
    AMREX_ROCSPARSE_SAFE_CALL
        (rocsparse_create_csr_descr(&r.mat_B, B.nrows, ncols, B.nnz,
                                    (void*)B.row_offset, (void*)B.col_index, (void*)B.mat,
                                    index_type, index_type, index_base, data_type));
    AMREX_ROCSPARSE_SAFE_CALL
        (rocsparse_create_csr_descr(&r.mat_C, nrows, ncols, 0,
                                    (void*)C.row_offset.data(), nullptr, nullptr,
                                    index_type, index_type, index_base, data_type));
    // D is unused because beta is zero, but a valid descriptor is required.
    AMREX_ROCSPARSE_SAFE_CALL
        (rocsparse_create_csr_descr(&r.mat_D, 0, 0, 0, nullptr, nullptr, nullptr,
                                    index_type, index_type, index_base, data_type));

    T alpha = T(1);
    T beta  = T(0);
    auto const op = rocsparse_operation_none;
    auto const alg = rocsparse_spgemm_alg_default;

    std::size_t buffer_size = 0;
    spgemm_rocsparse_check
        (rocsparse_spgemm(r.handle, op, op, &alpha, r.mat_A, r.mat_B, &beta, r.mat_D, r.mat_C,
                          data_type, alg, rocsparse_spgemm_stage_buffer_size,
                          &buffer_size, nullptr));
    r.buffer = (void*)The_Arena()->alloc(buffer_size);

    // This stage fills C's row offsets.
    spgemm_rocsparse_check
        (rocsparse_spgemm(r.handle, op, op, &alpha, r.mat_A, r.mat_B, &beta, r.mat_D, r.mat_C,
                          data_type, alg, rocsparse_spgemm_stage_nnz,
                          &buffer_size, r.buffer));

    std::int64_t c_nrows, c_ncols, c_nnz;
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_spmat_get_size(r.mat_C, &c_nrows, &c_ncols, &c_nnz));
    AMREX_ALWAYS_ASSERT(c_nrows == nrows && c_ncols == ncols);
    AMREX_ALWAYS_ASSERT(c_nnz < std::int64_t(std::numeric_limits<int>::max()));

    C.mat.resize(c_nnz);
    C.col_index.resize(c_nnz);
    C.nnz = c_nnz;
    AMREX_ROCSPARSE_SAFE_CALL
        (rocsparse_csr_set_pointers(r.mat_C, (void*)C.row_offset.data(),
                                    (void*)C.col_index.data(), (void*)C.mat.data()));

    spgemm_rocsparse_check
        (rocsparse_spgemm(r.handle, op, op, &alpha, r.mat_A, r.mat_B, &beta, r.mat_D, r.mat_C,
                          data_type, alg, rocsparse_spgemm_stage_compute,
                          &buffer_size, r.buffer));

    Gpu::streamSynchronize();
    The_Arena()->free(r.buffer);
    r.buffer = nullptr;

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

    // Releases everything, also when an OutOfMemoryError is thrown.
    struct Resources {
        sycl::queue& q;
        mkl::sparse::matrix_handle_t hA{}, hB{}, hC{};
        mkl::sparse::matmat_descr_t descr = nullptr;
        std::int64_t* size_buf = nullptr;
        void* buffer1 = nullptr;
        void* buffer2 = nullptr;
        explicit Resources (sycl::queue& a_q) : q(a_q) {}
        Resources (Resources const&) = delete;
        Resources& operator= (Resources const&) = delete;
        ~Resources () {
            q.wait();
            if (descr) { mkl::sparse::release_matmat_descr(&descr); }
            mkl::sparse::release_matrix_handle(q, &hA);
            mkl::sparse::release_matrix_handle(q, &hB);
            mkl::sparse::release_matrix_handle(q, &hC).wait();
            Gpu::streamSynchronize();
            if (buffer1) { The_Arena()->free(buffer1); }
            if (buffer2) { The_Arena()->free(buffer2); }
            if (size_buf) { The_Pinned_Arena()->free(size_buf); }
        }
    } r(q);

    mkl::sparse::init_matrix_handle(&r.hA);
    mkl::sparse::init_matrix_handle(&r.hB);
    mkl::sparse::init_matrix_handle(&r.hC);

#if defined(INTEL_MKL_VERSION) && (INTEL_MKL_VERSION < 20250300)
    mkl::sparse::set_csr_data(q, r.hA, I(nrows), I(B.nrows), mkl::index_base::zero,
                              (I*)A.row_offset, (I*)A.col_index, (T*)A.mat);
    mkl::sparse::set_csr_data(q, r.hB, I(B.nrows), I(ncols), mkl::index_base::zero,
                              (I*)B.row_offset, (I*)B.col_index, (T*)B.mat);
    mkl::sparse::set_csr_data(q, r.hC, I(nrows), I(ncols), mkl::index_base::zero,
                              C.row_offset.data(), dummy_col.data(), dummy_mat.data());
#else
    mkl::sparse::set_csr_data(q, r.hA, I(nrows), I(B.nrows), I(A.nnz), mkl::index_base::zero,
                              (I*)A.row_offset, (I*)A.col_index, (T*)A.mat);
    mkl::sparse::set_csr_data(q, r.hB, I(B.nrows), I(ncols), I(B.nnz), mkl::index_base::zero,
                              (I*)B.row_offset, (I*)B.col_index, (T*)B.mat);
    mkl::sparse::set_csr_data(q, r.hC, I(nrows), I(ncols), I(0), mkl::index_base::zero,
                              C.row_offset.data(), dummy_col.data(), dummy_mat.data());
#endif

    mkl::sparse::init_matmat_descr(&r.descr);
    mkl::sparse::set_matmat_data(r.descr,
                                 mkl::sparse::matrix_view_descr::general,
                                 mkl::transpose::nontrans,
                                 mkl::sparse::matrix_view_descr::general,
                                 mkl::transpose::nontrans,
                                 mkl::sparse::matrix_view_descr::general);

    using req = mkl::sparse::matmat_request;
    r.size_buf = (std::int64_t*)The_Pinned_Arena()->alloc(sizeof(std::int64_t));
    auto* size_buf = r.size_buf;

    mkl::sparse::matmat(q, r.hA, r.hB, r.hC, req::get_work_estimation_buf_size, r.descr,
                        size_buf, nullptr, {}).wait();
    r.buffer1 = (void*)The_Arena()->alloc(std::size_t(*size_buf));
    mkl::sparse::matmat(q, r.hA, r.hB, r.hC, req::work_estimation, r.descr,
                        size_buf, r.buffer1, {}).wait();

    mkl::sparse::matmat(q, r.hA, r.hB, r.hC, req::get_compute_buf_size, r.descr,
                        size_buf, nullptr, {}).wait();
    r.buffer2 = (void*)The_Arena()->alloc(std::size_t(*size_buf));
    mkl::sparse::matmat(q, r.hA, r.hB, r.hC, req::compute, r.descr,
                        size_buf, r.buffer2, {}).wait();

    mkl::sparse::matmat(q, r.hA, r.hB, r.hC, req::get_nnz, r.descr,
                        size_buf, nullptr, {}).wait();
    Long const c_nnz = *size_buf;
    AMREX_ALWAYS_ASSERT(c_nnz < Long(std::numeric_limits<I>::max()));

    C.mat.resize(c_nnz);
    C.col_index.resize(c_nnz);
    C.nnz = c_nnz;
#if defined(INTEL_MKL_VERSION) && (INTEL_MKL_VERSION < 20250300)
    mkl::sparse::set_csr_data(q, r.hC, I(nrows), I(ncols), mkl::index_base::zero,
                              C.row_offset.data(), C.col_index.data(), C.mat.data());
#else
    mkl::sparse::set_csr_data(q, r.hC, I(nrows), I(ncols), I(c_nnz), mkl::index_base::zero,
                              C.row_offset.data(), C.col_index.data(), C.mat.data());
#endif

    mkl::sparse::matmat(q, r.hA, r.hB, r.hC, req::finalize, r.descr,
                        size_buf, nullptr, {}).wait();

    Gpu::streamSynchronize();
    The_Arena()->free(r.buffer1);
    The_Arena()->free(r.buffer2);
    r.buffer1 = nullptr;
    r.buffer2 = nullptr;

    C.sort(); // oneMKL does not sort the output.

    return C;
}

#endif

#ifdef AMREX_USE_GPU

/**
 * Low-memory fallback: `f` multiplies blocks of A's rows by B one at a time,
 * and a block that runs out of memory is split in half.
 */
template <typename T, template<typename> class V, typename I>
CSR<T,V,I> spgemm_local_chunked (Long nrows, Long ncols,
                                 CsrView<T const,I> const& A, CsrView<T const,I> const& B,
                                 CSR<T,V,I> (*f) (Long, Long, CsrView<T const,I> const&,
                                                  CsrView<T const,I> const&))
{
    struct Chunk {
        Long r0;
        Long r1;
        CSR<T,V,I> C;
    };
    Vector<Chunk> done;

    // Row ranges to do, processed from the back so that `done` is in row order.
    Vector<std::pair<Long,Long>> todo;
    todo.emplace_back(nrows/2, nrows);
    todo.emplace_back(0, nrows/2);

    while (!todo.empty()) {
        auto const [r0, r1] = todo.back();
        todo.pop_back();
        Long const n = r1 - r0;
        if (n == 0) { continue; }

        I a[2];
        Gpu::dtoh_memcpy_async(a  , A.row_offset+r0, sizeof(I));
        Gpu::dtoh_memcpy_async(a+1, A.row_offset+r1, sizeof(I));
        Gpu::streamSynchronize();
        I const a0 = a[0];
        Long const annz = Long(a[1]) - Long(a0);

        try {
            if (annz == 0) {
                done.push_back(Chunk{r0, r1, spgemm_empty<T,V,I>(n)});
            } else {
                V<I> row(n+1);
                auto* prow = row.data();
                auto const* arow = A.row_offset + r0;
                ParallelFor(n+1, [=] AMREX_GPU_DEVICE (Long i) { prow[i] = arow[i] - a0; });
                CsrView<T const,I> Ab{A.mat + a0, A.col_index + a0, row.data(), annz, n};
                done.push_back(Chunk{r0, r1, f(n, ncols, Ab, B)});
            }
        } catch (OutOfMemoryError const&) {
            if (n == 1) { throw; }
            todo.emplace_back(r0 + n/2, r1);
            todo.emplace_back(r0, r0 + n/2);
        }
    }

    Long total = 0;
    for (auto const& c : done) { total += c.C.nnz; }
    AMREX_ALWAYS_ASSERT(total <= Long(std::numeric_limits<I>::max()));

    CSR<T,V,I> C;
    C.resize(nrows, total);
    Long p = 0;
    for (auto& c : done) {
        auto* crow = C.row_offset.data() + c.r0;
        auto const* src = c.C.row_offset.data();
        auto const offset = I(p);
        ParallelFor(c.r1 - c.r0 + 1, [=] AMREX_GPU_DEVICE (Long i) {
            crow[i] = src[i] + offset;
        });
        Gpu::copyAsync(Gpu::deviceToDevice, c.C.col_index.begin(), c.C.col_index.end(),
                       C.col_index.begin() + p);
        Gpu::copyAsync(Gpu::deviceToDevice, c.C.mat.begin(), c.C.mat.end(),
                       C.mat.begin() + p);
        p += c.C.nnz;
        c.C = CSR<T,V,I>{}; // stream-ordered free
    }
    Gpu::streamSynchronize();

    return C;
}

#endif

/**
 * \brief Local CSR product C = A * B.
 *
 * `A` must have `nrows` rows and column indices in `[0, B.nrows)`; `B` must
 * have column indices in `[0, ncols)`. Rows of the result are sorted. On
 * GPUs, a slower algorithm that uses less memory is used if the default one
 * runs out of memory.
 */
template <typename T, template<typename> class V, typename I>
CSR<T,V,I> spgemm_local (Long nrows, Long ncols,
                         CsrView<T const,I> const& A, CsrView<T const,I> const& B)
{
    AMREX_ASSERT(A.nrows == nrows);

    if (nrows <= 0 || ncols <= 0 || A.nnz <= 0 || B.nnz <= 0 || B.nrows <= 0) {
        return spgemm_empty<T,V,I>(nrows);
    }

#if !defined(AMREX_USE_GPU)
    return spgemm_local_cpu<T,V,I>(nrows, ncols, A, B);
#elif defined(AMREX_USE_CUDA)
    try {
        return spgemm_local_cusparse<T,V,I>(nrows, ncols, A, B);
    } catch (OutOfMemoryError const&) {
        return spgemm_local_chunked<T,V,I>(nrows, ncols, A, B,
                                           spgemm_local_cusparse<T,V,I>);
    }
#elif defined(AMREX_USE_HIP)
    try {
        return spgemm_local_rocsparse<T,V,I>(nrows, ncols, A, B);
    } catch (OutOfMemoryError const&) {
        return spgemm_local_chunked<T,V,I>(nrows, ncols, A, B,
                                           spgemm_local_rocsparse<T,V,I>);
    }
#elif defined(AMREX_USE_SYCL)
    try {
        return spgemm_local_onemkl<T,V,I>(nrows, ncols, A, B);
    } catch (OutOfMemoryError const&) {
        return spgemm_local_chunked<T,V,I>(nrows, ncols, A, B,
                                           spgemm_local_onemkl<T,V,I>);
    }
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
template <typename T, template <typename> class Allocator>
SpMatrix<T,Allocator>
SpGEMM (SpMatrix<T,Allocator> const& A, SpMatrix<T,Allocator> const& B,
        AlgPartition const& col_partition)
{
    BL_PROFILE("SpGEMM");

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
        (nrows, ncols, Am.m_csr_local.const_view(), Bm.m_csr_local.const_view());
    SpMat C;
    C.define_split(Am.partition(), col_partition, std::move(Ch), ncols, nullptr, 0);
    return C;

#endif
}

}

#endif
