#ifndef AMREX_SPMV_H_
#define AMREX_SPMV_H_
#include <AMReX_Config.H>

#include <AMReX_AlgVector.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_SpMatrix.H>

#if defined(AMREX_USE_CUDA)
#  include <cusparse.h>
#elif defined(AMREX_USE_HIP)
#  include <hipsparse/hipsparse.h>
#elif defined(AMREX_USE_DPCPP)
#  include <oneapi/mkl/spblas.hpp>
#endif

namespace amrex {

template <typename T>
void SpMV (AlgVector<T>& y, SpMatrix<T> const& A, AlgVector<T> const& x)
{
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

    hipsparseHandle_t handle;
    hipsparseCreate(&handle);
    hipsparseSetStream(handle, Gpu::gpuStream());

    hipDataType data_type;
    if constexpr (std::is_same_v<T,float>) {
        data_type = HIP_R_32F;
    } else if constexpr (std::is_same_v<T,double>) {
        data_type = HIP_R_64F;
    } else if constexpr (std::is_same_v<T,GpuComplex<float>>) {
        data_type = HIP_C_32F;
    } else if constexpr (std::is_same_v<T,GpuComplex<double>>) {
        data_type = HIP_C_64F;
    } else {
        amrex::Abort("SpMV: unsupported data type");
    }

    hipsparseIndexType_t index_type =  HIPSPARSE_INDEX_64I;

    hipsparseSpMatDescr_t mat_descr;
    hipsparseCreateCsr(&mat_descr, nrows, ncols, nnz, (void*)row, (void*)col, (void*)mat,
                       index_type, index_type, HIPSPARSE_INDEX_BASE_ZERO, data_type);

    hipsparseDnVecDescr_t x_descr;
    hipsparseCreateDnVec(&x_descr, ncols, (void*)px, data_type);

    hipsparseDnVecDescr_t y_descr;
    hipsparseCreateDnVec(&y_descr, nrows, (void*)py, data_type);

    T alpha = T(1.0);
    T beta = T(0.0);

    std::size_t buffer_size;
    hipsparseSpMV_bufferSize(handle, HIPSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat_descr, x_descr,
                             &beta, y_descr, data_type, HIPSPARSE_SPMV_ALG_DEFAULT, &buffer_size);

    void* pbuffer = (void*)The_Arena()->alloc(buffer_size);

    hipsparseSpMV(handle, HIPSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat_descr, x_descr,
                  &beta, y_descr, data_type, HIPSPARSE_SPMV_ALG_DEFAULT, pbuffer);

    Gpu::streamSynchronize();

    hipsparseDestroySpMat(mat_descr);
    hipsparseDestroyDnVec(x_descr);
    hipsparseDestroyDnVec(y_descr);
    hipsparseDestroy(handle);
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
