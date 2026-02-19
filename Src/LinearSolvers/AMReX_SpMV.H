#ifndef AMREX_SPMV_H_
#define AMREX_SPMV_H_
#include <AMReX_Config.H>

#include <AMReX_AlgVector.H>
#include <AMReX_AlgVecUtil.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_SpMatrix.H>

namespace amrex {

template <typename T>
void SpMV (Long nrows, Long ncols, T* AMREX_RESTRICT py, CsrView<T const> const& A,
           T const* AMREX_RESTRICT px)
{
    T const* AMREX_RESTRICT mat = A.mat;
    auto const* AMREX_RESTRICT col = A.col_index;
    auto const* AMREX_RESTRICT row = A.row_offset;

#if defined(AMREX_USE_GPU)

    Long const nnz = A.nnz;

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
        amrex::Abort("SpMV: unsupported data type");
    }

    cusparseIndexType_t index_type = CUSPARSE_INDEX_64I;

    cusparseSpMatDescr_t mat_descr;
    AMREX_CUSPARSE_SAFE_CALL
        (cusparseCreateCsr(&mat_descr, nrows, ncols, nnz,
                           (void*)row, (void*)col, (void*)mat,
                           index_type, index_type, CUSPARSE_INDEX_BASE_ZERO,
                           data_type));

    cusparseDnVecDescr_t x_descr;
    AMREX_CUSPARSE_SAFE_CALL(cusparseCreateDnVec(&x_descr, ncols, (void*)px, data_type));

    cusparseDnVecDescr_t y_descr;
    AMREX_CUSPARSE_SAFE_CALL(cusparseCreateDnVec(&y_descr, nrows, (void*)py, data_type));

    T alpha = T(1);
    T beta = T(0);

    std::size_t buffer_size;
    AMREX_CUSPARSE_SAFE_CALL
        (cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha, mat_descr, x_descr, &beta, y_descr,
                                 data_type, CUSPARSE_SPMV_ALG_DEFAULT,
                                 &buffer_size));

    auto* pbuffer = (void*)The_Arena()->alloc(buffer_size);

    AMREX_CUSPARSE_SAFE_CALL
        (cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                      &alpha, mat_descr, x_descr, &beta, y_descr,
                      data_type, CUSPARSE_SPMV_ALG_DEFAULT, pbuffer));

    Gpu::streamSynchronize();

    AMREX_CUSPARSE_SAFE_CALL(cusparseDestroySpMat(mat_descr));
    AMREX_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(x_descr));
    AMREX_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(y_descr));
    AMREX_CUSPARSE_SAFE_CALL(cusparseDestroy(handle));
    The_Arena()->free(pbuffer);

#elif defined(AMREX_USE_HIP)

    rocsparse_handle handle;
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_create_handle(&handle));
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_set_stream(handle, Gpu::gpuStream()));

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
        amrex::Abort("SpMV: unsupported data type");
    }

    rocsparse_indextype index_type = rocsparse_indextype_i64;

    rocsparse_spmat_descr mat_descr;
    AMREX_ROCSPARSE_SAFE_CALL(
        rocsparse_create_csr_descr(&mat_descr, nrows, ncols, nnz,
                                   (void*)row, (void*)col, (void*)mat,
                                   index_type, index_type,
                                   rocsparse_index_base_zero, data_type));

    rocsparse_dnvec_descr x_descr;
    AMREX_ROCSPARSE_SAFE_CALL(
        rocsparse_create_dnvec_descr(&x_descr, ncols, (void*)px, data_type));

    rocsparse_dnvec_descr y_descr;
    AMREX_ROCSPARSE_SAFE_CALL(
        rocsparse_create_dnvec_descr(&y_descr, nrows, (void*)py, data_type));

    T alpha = T(1.0);
    T beta = T(0.0);

#if (HIP_VERSION_MAJOR >= 7)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

    std::size_t buffer_size;
    auto err0 =
        rocsparse_spmv(handle, rocsparse_operation_none, &alpha, mat_descr, x_descr,
                       &beta, y_descr, data_type, rocsparse_spmv_alg_default,
#if (HIP_VERSION_MAJOR >= 6)
                       rocsparse_spmv_stage_buffer_size,
#endif
                       &buffer_size, nullptr);
    AMREX_ROCSPARSE_SAFE_CALL(err0);

    void* pbuffer = (void*)The_Arena()->alloc(buffer_size);

#if (HIP_VERSION_MAJOR >= 6)
    AMREX_ROCSPARSE_SAFE_CALL(
        rocsparse_spmv(handle, rocsparse_operation_none, &alpha, mat_descr, x_descr,
                       &beta, y_descr, data_type, rocsparse_spmv_alg_default,
                       rocsparse_spmv_stage_preprocess, &buffer_size, pbuffer));
#endif

    auto err1 =
        rocsparse_spmv(handle, rocsparse_operation_none, &alpha, mat_descr, x_descr,
                       &beta, y_descr, data_type, rocsparse_spmv_alg_default,
#if (HIP_VERSION_MAJOR >= 6)
                       rocsparse_spmv_stage_compute,
#endif
                       &buffer_size, pbuffer);
    AMREX_ROCSPARSE_SAFE_CALL(err1);

#if (HIP_VERSION_MAJOR >= 7)
// xxxxx HIP TODO: rocsparse_spmv has been deprecated.
#pragma clang diagnostic pop
#endif

    Gpu::streamSynchronize();

    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_destroy_spmat_descr(mat_descr));
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_destroy_dnvec_descr(x_descr));
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_destroy_dnvec_descr(y_descr));
    AMREX_ROCSPARSE_SAFE_CALL(rocsparse_destroy_handle(handle));
    The_Arena()->free(pbuffer);

#elif defined(AMREX_USE_SYCL)

    mkl::sparse::matrix_handle_t handle{};
    mkl::sparse::init_matrix_handle(&handle);

#if defined(INTEL_MKL_VERSION) && (INTEL_MKL_VERSION < 20250300)
    amrex::ignore_unused(nnz);
    mkl::sparse::set_csr_data(Gpu::Device::streamQueue(), handle, nrows, ncols,
                              mkl::index_base::zero, (Long*)row, (Long*)col, (T*)mat);
#else
    mkl::sparse::set_csr_data(Gpu::Device::streamQueue(), handle, nrows, ncols, nnz,
                              mkl::index_base::zero, (Long*)row, (Long*)col, (T*)mat);
#endif
    mkl::sparse::gemv(Gpu::Device::streamQueue(), mkl::transpose::nontrans,
                      T(1), handle, px, T(0), py);

    auto ev = mkl::sparse::release_matrix_handle(Gpu::Device::streamQueue(), &handle);
    ev.wait();

#endif

    AMREX_GPU_ERROR_CHECK();

#else

    amrex::ignore_unused(ncols);

#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (Long i = 0; i < nrows; ++i) {
        T r = 0;
        for (Long j = row[i]; j < row[i+1]; ++j) {
            r += mat[j] * px[col[j]];
        }
        py[i] = r;
    }

#endif
}

template <typename T, template<typename> class AllocM, typename AllocV>
void SpMV (AlgVector<T,AllocV>& y, SpMatrix<T,AllocM> const& A,
           AlgVector<T,AllocV> const& x)
{
    // xxxxx TODO: Is it worth to cache cusparse and rocsparse handles?

    const_cast<SpMatrix<T,AllocM>&>(A).startComm_mv(x);

    // Diagonal part
    SpMV<T>(y.numLocalRows(), x.numLocalRows(), y.data(), A.m_csr.const_view(),
            x.data());

    const_cast<SpMatrix<T,AllocM>&>(A).finishComm_mv(y);
}

//! res = b - A*x
template <typename T, template<typename> class AllocM, typename AllocV>
void computeResidual (AlgVector<T,AllocV>& res, SpMatrix<T,AllocM> const& A,
                      AlgVector<T,AllocV> const& x, AlgVector<T,AllocV> const& b)
{
    SpMV(res, A, x);
    Xpay(res, T(-1), b);
}

}

#endif
