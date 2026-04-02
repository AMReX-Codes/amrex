#ifndef AMREX_SPMAT_UTIL_H_
#define AMREX_SPMAT_UTIL_H_
#include <AMReX_Config.H>

#include <AMReX_SpMatrix.H>

namespace amrex {

/**
 * \file AMReX_SpMatUtil.H
 *
 * Utility helpers for sparse matrices (CSR transpose, etc.).
 */

/**
 * \brief Build the transpose CSR of \p csr.
 *
 * \param csr   Source matrix in CSR form (local portion).
 * \param ncols Total number of columns in the source matrix.
 * \return CSR representing the transpose.
 */
template <typename T, template <typename> class V>
CSR<T,V> transpose (CSR<T,V> const& csr, Long ncols)
{
    Long nnz = csr.nnz;

    CSR<T,V> csrt;
    csrt.mat.resize(nnz);
    csrt.col_index.resize(nnz);
    csrt.row_offset.resize(ncols+1);
    csrt.nnz = nnz;

    detail::transpose(csrt.view(), csr.const_view());

    return csrt;
}

/**
 * \brief Transpose a distributed SpMatrix, producing a new matrix with column partition \p col_partition.
 *
 * \param A             Input matrix.
 * \param col_partition Desired column partition (becomes rows of the transpose).
 * \return Transposed sparse matrix.
 */
template <typename T, template<typename> class Allocator>
SpMatrix<T,Allocator> transpose (SpMatrix<T,Allocator> const& A,
                                 AlgPartition col_partition)
{
    const_cast<SpMatrix<T,Allocator>&>(A).startComm_tr(col_partition);

    SpMatrix<T,Allocator> AT;
    AT.m_row_begin = col_partition[ParallelDescriptor::MyProc()];
    AT.m_row_end   = col_partition[ParallelDescriptor::MyProc()+1];
    AT.m_partition = std::move(col_partition);
    AT.m_nnz       = A.m_csr.nnz;

    // transpose the diagonal part
    auto at_csr = transpose(A.m_csr, AT.numLocalRows());
    AT.m_csr = std::move(at_csr);

    const_cast<SpMatrix<T,Allocator>&>(A).finishComm_tr(AT);

    return AT;
}

}

#endif
