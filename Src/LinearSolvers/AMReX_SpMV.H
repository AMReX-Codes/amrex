#ifndef AMREX_SPMV_H_
#define AMREX_SPMV_H_
#include <AMReX_Config.H>

#include <AMReX_AlgVector.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_SpMatrix.H>

#if defined(AMREX_USE_CUDA)
#  include <cusparse.h>
#elif defined(AMREX_USE_HIP)
#  include <rocsparse/rocsparse.h>
#elif defined(AMREX_USE_DPCPP)
#  include <oneapi/mkl/spblas.hpp>
#endif

namespace amrex {

template <typename T>
void SpMV (AlgVector<T>& y, SpMatrix<T> const& A, AlgVector<T> const& x)
{
    // xxxxx TODOL We might want to cache the cusparse and rocsparse handles

    // xxxxx TODO: let's assume it's square matrix for now.
    AMREX_ALWAYS_ASSERT(x.partition() == y.partition() &&
                        x.partition() == A.partition());

    const_cast<SpMatrix<T>&>(A).startComm(x);

    T      * AMREX_RESTRICT py = y.data();
    T const* AMREX_RESTRICT px = x.data();
    T const* AMREX_RESTRICT mat = A.data();
    auto const* AMREX_RESTRICT col = A.columnIndex();
    auto const* AMREX_RESTRICT row = A.rowOffset();

#if defined(AMREX_USE_GPU)

    Long const nrows = A.numLocalRows();
    Long const ncols = x.numLocalRows();
    Long const nnz = A.numLocalNonZero();

#if defined(AMREX_USE_CUDA)

    cusparseHandle_t handle;
    cusparseCreate(&handle);
    cusparseSetStream(handle, Gpu::gpuStream());

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
    cusparseCreateCsr(&mat_descr, nrows, ncols, nnz, (void*)row, (void*)col, (void*)mat,
                      index_type, index_type, CUSPARSE_INDEX_BASE_ZERO, data_type);

    cusparseDnVecDescr_t x_descr;
    cusparseCreateDnVec(&x_descr, ncols, (void*)px, data_type);

    cusparseDnVecDescr_t y_descr;
    cusparseCreateDnVec(&y_descr, nrows, (void*)py, data_type);

    T alpha = T(1);
    T beta = T(0);

    std::size_t buffer_size;
    cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat_descr, x_descr,
                            &beta, y_descr, data_type, CUSPARSE_SPMV_ALG_DEFAULT, &buffer_size);

    auto* pbuffer = (void*)The_Arena()->alloc(buffer_size);

    cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat_descr, x_descr,
                 &beta, y_descr, data_type, CUSPARSE_SPMV_ALG_DEFAULT, pbuffer);

    Gpu::streamSynchronize();

    cusparseDestroySpMat(mat_descr);
    cusparseDestroyDnVec(x_descr);
    cusparseDestroyDnVec(y_descr);
    cusparseDestroy(handle);
    The_Arena()->free(pbuffer);

#elif defined(AMREX_USE_HIP)

    rocsparse_handle handle;
    rocsparse_create_handle(&handle);
    rocsparse_set_stream(handle, Gpu::gpuStream());

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
    rocsparse_create_csr_descr(&mat_descr, nrows, ncols, nnz, (void*)row, (void*)col,
                               (void*)mat, index_type, index_type,
                               rocsparse_index_base_zero, data_type);

    rocsparse_dnvec_descr x_descr;
    rocsparse_create_dnvec_descr(&x_descr, ncols, (void*)px, data_type);

    rocsparse_dnvec_descr y_descr;
    rocsparse_create_dnvec_descr(&y_descr, nrows, (void*)py, data_type);

    T alpha = T(1.0);
    T beta = T(0.0);

    std::size_t buffer_size;
    rocsparse_spmv(handle, rocsparse_operation_none, &alpha, mat_descr, x_descr,
                   &beta, y_descr, data_type, rocsparse_spmv_alg_default,
#if (HIP_VERSION_MAJOR >= 6)
                   rocsparse_spmv_stage_buffer_size,
#endif
                   &buffer_size, nullptr);

    void* pbuffer = (void*)The_Arena()->alloc(buffer_size);

#if (HIP_VERSION_MAJOR >= 6)
    rocsparse_spmv(handle, rocsparse_operation_none, &alpha, mat_descr, x_descr,
                   &beta, y_descr, data_type, rocsparse_spmv_alg_default,
                   rocsparse_spmv_stage_preprocess, &buffer_size, pbuffer);
#endif

    rocsparse_spmv(handle, rocsparse_operation_none, &alpha, mat_descr, x_descr,
                   &beta, y_descr, data_type, rocsparse_spmv_alg_default,
#if (HIP_VERSION_MAJOR >= 6)
                   rocsparse_spmv_stage_compute,
#endif
                   &buffer_size, pbuffer);

    Gpu::streamSynchronize();

    rocsparse_destroy_spmat_descr(mat_descr);
    rocsparse_destroy_dnvec_descr(x_descr);
    rocsparse_destroy_dnvec_descr(y_descr);
    rocsparse_destroy_handle(handle);
    The_Arena()->free(pbuffer);

#elif defined(AMREX_USE_DPCPP)

    mkl::sparse::matrix_handle_t handle{};
    mkl::sparse::set_csr_data(Gpu::Device::streamQueue(), handle, nrows, ncols, mkl::index_base::zero,
                              (Long*)row, (Long*)col, (T*)mat);
    mkl::sparse::gemv(Gpu::Device::streamQueue(), mkl::transpose::nontrans,
                      T(1), handle, px, T(0), py);

#endif

    AMREX_GPU_ERROR_CHECK();

#else

    Long const ny = y.numLocalRows();
    for (Long i = 0; i < ny; ++i) {
        T r = 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel for reduction(+:r)
#endif
        for (Long j = row[i]; j < row[i+1]; ++j) {
            r += mat[j] * px[col[j]];
        }
        py[i] = r;
    }

#endif

    const_cast<SpMatrix<T>&>(A).finishComm(y);
}

}

#endif
