#ifndef AMREX_SPMAT_UTIL_H_
#define AMREX_SPMAT_UTIL_H_
#include <AMReX_Config.H>

#include <AMReX_BLProfiler.H>
#include <AMReX_SpMatrix.H>
#include <AMReX_Random.H>

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
template <typename IO = Long, typename T, template <typename> class V, typename I>
CSR<T,V,IO> transpose (CSR<T,V,I> const& csr, Long ncols)
{
    Long nnz = csr.nnz;

    CSR<T,V,IO> csrt;
    csrt.mat.resize(nnz);
    csrt.col_index.resize(nnz);
    csrt.row_offset.resize(ncols+1);
    csrt.nnz = nnz;

    detail::transpose(csrt.view(), csr.const_view());

    return csrt;
}

/**
 * \brief Transpose sparse matrix
 *
 * \param A Sparse matrix
 * \param col_partition Column partition of `A` (row partition of `A^T`).
 *
 * \return Transposed matrix
 */
template <typename T, template<typename> class Allocator>
SpMatrix<T,Allocator> transpose (SpMatrix<T,Allocator> const& A,
                                 AlgPartition const& col_partition)
{
    BL_PROFILE("SpMatrix::transpose");

    auto& Am = const_cast<SpMatrix<T,Allocator>&>(A);
    Am.setColumnPartition(col_partition);
    Am.startComm_tr(col_partition);

    // The transpose is built directly in split form: its local block is
    // the transpose of A's local block (columns are A's local row
    // indices), and its remote block comes from the other ranks.
    SpMatrix<T,Allocator> AT;
    AT.m_row_begin = col_partition[ParallelDescriptor::MyProc()];
    AT.m_row_end   = col_partition[ParallelDescriptor::MyProc()+1];
    AT.m_partition = col_partition;
    AT.m_nnz       = A.m_csr_local.nnz;
    AT.m_csr_local = transpose<int>(A.m_csr_local, AT.numLocalRows());
    AT.m_col_partition = A.partition();
    AT.m_col_begin = A.globalRowBegin();
    AT.m_col_end   = A.globalRowEnd();
    AT.m_split = true;
#ifdef AMREX_USE_MPI
    AT.m_ri_ltor.resize(AT.numLocalRows(), -1);
    AT.m_remote_cols_vv.resize(ParallelDescriptor::NProcs());
#endif

    Am.finishComm_tr(AT);

    return AT;
}

//! Return an identity matrix with given partition
template <typename T, template<typename> class Allocator = DefaultAllocator>
SpMatrix<T,Allocator> IdentityMatrix (AlgPartition partition)
{
    SpMatrix<T,Allocator> I(std::move(partition), 1);
    I.setVal([=] AMREX_GPU_DEVICE (Long row, Long* col, T* val) {
        *col = row; *val = T(1);
    }, CsrSorted{true});
    return I;
}

/**
 * \brief Return a random sparse matrix
 *
 * This function generates a random sparse matrix. For each row, the number
 * of nonzeros is drawn from a Poisson distribution with parameter `lambda`
 * and capped by `nnz_per_row_max`. Column indices are chosen randomly.
 * Matrix values are random numbers uniformly distributed in the range of
 * [0,1).
 *
 * \param partition Parallel partition information.
 * \param nrows Total number of rows.
 * \param ncols Total number of columns.
 * \param lambda Parameter for the Poisson distribution.
 * \param nnz_per_row_max Maximum number of nonzeros in a row.
 */
template <typename T, template<typename> class Allocator = DefaultAllocator>
SpMatrix<T,Allocator> RandomMatrix (AlgPartition partition, Long nrows, Long ncols,
                                    Real lambda, int nnz_per_row_max)
{
    AMREX_ALWAYS_ASSERT(partition.numGlobalRows() == nrows);

    auto nlrows = partition.numLocalRows();
    Long lnnz = nlrows * nnz_per_row_max;
    typename SpMatrix<T,Allocator>::template container_type<T> mat(lnnz);
    typename SpMatrix<T,Allocator>::template container_type<Long> col_index(lnnz);
    typename SpMatrix<T,Allocator>::template container_type<Long> row_offset(nlrows+1, 0);
    auto* pmat = mat.data();
    auto* pcol = col_index.data();
    auto* prow = row_offset.data();

    auto ncols_max = (unsigned int)
        std::min(ncols, Long(std::numeric_limits<unsigned int>::max()));
    // A row cannot have more distinct columns than ncols.
    auto nnz_row_cap = int(std::min(Long(nnz_per_row_max), ncols));

    ParallelForRNG(nlrows, [=] AMREX_GPU_DEVICE (Long lrow, RandomEngine const& eng)
    {
        auto os = lrow * nnz_per_row_max;
        prow[lrow] = os;
        if (lrow+1 == nlrows) {
            prow[lrow+1] = lnnz;
        }

        auto nnz_row = int(std::min(Long(amrex::RandomPoisson(lambda, eng)),
                                    Long(nnz_row_cap)));

        for (int j = 0; j < nnz_row; ++j) {
            pmat[os+j] = Random(eng);
            while (true) {
                auto c = Long(Random_int(ncols_max, eng));
                bool new_c = true;
                for (int jj = 0; jj < j; ++jj) {
                    new_c = new_c && (c != pcol[os+jj]);
                }
                if (new_c) {
                    pcol[os+j] = c;
                    break;
                }
            }
        }

        for (int j = nnz_row; j < nnz_per_row_max; ++j) {
            pcol[os+j] = -1;
            pmat[os+j] = T(0);
        }
    });

    SpMatrix<T,Allocator> r;
    r.define(std::move(partition), pmat, pcol, lnnz, prow,
             CsrSorted{false}, CsrValid{false});

    return r;
}

/**
 * \brief Check whether two sparse matrices are almost equal
 *
 * This function returns true if and only if the two matrices have the same
 * partition and sparse pattern, and all corresponding matrix entries are
 * almost equal within a specified tolerance measured in units in the last
 * place (ULPs).
 *
 * The comparison is performed locally. No MPI communication is
 * involved. Thus, the result may be different among processes; reduce it
 * with `ParallelDescriptor::ReduceBoolAnd` if a global answer is needed.
 *
 * If only one matrix has a column partition, the other one is given the
 * same column partition, which it then keeps. If neither has one, the
 * unsplit matrices are compared.
 *
 * \param A First matrix.
 * \param B Second matrix.
 * \param ulp Maximum allowed difference in ULPs.
 *
 * \return `true` if the matrices ae almost equal, `false` otherwise.
 */
template <typename T, template<typename> class Allocator = DefaultAllocator>
bool almostEqual (SpMatrix<T,Allocator> const& A, SpMatrix<T,Allocator> const& B,
                  int ulp = 2)
{
    auto const& col_partition_A = A.columnPartition();
    auto const& col_partition_B = B.columnPartition();
    if (col_partition_A.empty() && ! col_partition_B.empty()) {
        const_cast<SpMatrix<T,Allocator>&>(A).setColumnPartition(col_partition_B);
    } else if (! col_partition_A.empty() && col_partition_B.empty()) {
        const_cast<SpMatrix<T,Allocator>&>(B).setColumnPartition(col_partition_A);
    } else if (col_partition_A.empty() && col_partition_B.empty()) {
        // Neither has been split: compare the unsplit CSR with global columns.
        auto& a = const_cast<SpMatrix<T,Allocator>&>(A);
        auto& b = const_cast<SpMatrix<T,Allocator>&>(B);
        Long const nnz = a.numLocalNonZeros();
        Long const nrows = a.numLocalRows();
        if (a.globalRowBegin() != b.globalRowBegin() ||
            nnz != b.numLocalNonZeros() || nrows != b.numLocalRows())
        {
            return false;
        }
        T const* amat = a.data();
        T const* bmat = b.data();
        Long const* acol = a.columnIndex();
        Long const* bcol = b.columnIndex();
        Long const* arow = a.rowOffset();
        Long const* brow = b.rowOffset();
        Long N = std::max(nnz, nrows);
        auto r = Reduce::Min<int>(N, [=] AMREX_GPU_DEVICE (Long i) {
            bool t = true;
            if (i < nnz) {
                t = amrex::almostEqual(amat[i], bmat[i], ulp) && acol[i] == bcol[i];
            }
            if (t && i < nrows) {
                t = arow[i] == brow[i];
            }
            return int(t);
        }, 1);
        return r;
    }

    auto const& pca = A.const_parcsr();
    auto const& pcb = B.const_parcsr();

    if (pca.row_begin  != pcb.row_begin  ||
        pca.col_begin  != pcb.col_begin  ||
        pca.csr0.nnz   != pcb.csr0.nnz   ||
        pca.csr0.nrows != pcb.csr0.nrows ||
        pca.csr1.nnz   != pcb.csr1.nnz   ||
        pca.csr1.nrows != pcb.csr1.nrows)
    {
        return false;
    }

    Long N = std::max({pca.csr0.nnz, pca.csr0.nrows, pca.csr1.nnz, pca.csr1.nrows});
    auto r = Reduce::Min<int>(N, [=] AMREX_GPU_DEVICE (Long i) {
        bool t = true;
        if (i < pca.csr0.nnz) {
            t = amrex::almostEqual(pca.csr0.mat[i], pcb.csr0.mat[i], ulp) &&
                pca.csr0.col_index[i] == pcb.csr0.col_index[i];
        }
        if (t && i < pca.csr0.nrows) {
            t = pca.csr0.row_offset[i] == pcb.csr0.row_offset[i];
            // Which local rows have remote entries.
            if (t && pca.csr1.nnz > 0) {
                t = pca.row_map[i] == pcb.row_map[i];
            }
        }
        if (t && i < pca.csr1.nnz) {
            // Remote column indices are local; compare the global ones.
            t = amrex::almostEqual(pca.csr1.mat[i], pcb.csr1.mat[i], ulp) &&
                pca.col_map[pca.csr1.col_index[i]] == pcb.col_map[pcb.csr1.col_index[i]];
        }
        if (t && i < pca.csr1.nrows) {
            t = pca.csr1.row_offset[i] == pcb.csr1.row_offset[i];
        }
        return int(t);
    }, 1);

    return r;
}

}

#endif
