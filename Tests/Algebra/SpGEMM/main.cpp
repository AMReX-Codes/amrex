#include <AMReX_Algebra.H>
#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>

using namespace amrex;

namespace {

// Partition with uneven sizes; rank `empty_rank` (if valid) gets no rows.
AlgPartition make_partition (Long nrows, int empty_rank, int seed)
{
    int const nprocs = ParallelDescriptor::NProcs();
    Vector<Long> rows(nprocs+1, 0);
    Vector<Long> w(nprocs, 1);
    for (int i = 0; i < nprocs; ++i) {
        w[i] = 1 + (i*7 + seed) % 5;
        if (i == empty_rank && nprocs > 1) { w[i] = 0; }
    }
    Long wsum = std::accumulate(w.begin(), w.end(), Long(0)); // >= 1
    Long acc = 0;
    for (int i = 0; i < nprocs; ++i) {
        acc += w[i];
        rows[i+1] = nrows * acc / wsum;
    }
    return AlgPartition(rows);
}

bool all_true (bool b)
{
    ParallelDescriptor::ReduceBoolAnd(b);
    return b;
}

// Compare (A*B)x with A(Bx).
template <typename T>
bool check_spmv (SpMatrix<T> const& A, SpMatrix<T> const& B, SpMatrix<T> const& AB,
                 AlgPartition const& xpart, T tol)
{
    AlgVector<T> x(xpart);
    auto* px = x.data();
    Long begin = x.globalBegin();
    ParallelFor(x.numLocalRows(), [=] AMREX_GPU_DEVICE (Long i) {
        px[i] = T(1) + T(0.25) * T((i + begin) % 7);
    });

    AlgVector<T> t(B.partition());
    AlgVector<T> y1(A.partition());
    AlgVector<T> y2(A.partition());
    SpMV(t, B, x);
    SpMV(y1, A, t);
    SpMV(y2, AB, x);

    auto y1max = y1.norminf();
    auto* p1 = y1.data();
    auto const* p2 = y2.data();
    ParallelFor(y1.numLocalRows(), [=] AMREX_GPU_DEVICE (Long i) {
        p1[i] -= p2[i];
    });
    Gpu::streamSynchronize();
    auto err = y1.norminf();
    return err <= tol * y1max;
}

}

int main (int argc, char *argv[])
{
    amrex::Initialize(argc, argv);
    {
        int const nprocs = ParallelDescriptor::NProcs();

        // Identity Matrix
        for (int ipart = 0; ipart < 3; ++ipart)
        {
            Long nrows = 30;
            Long ncols = 27;
            AlgPartition rpart = (ipart == 0) ? AlgPartition(nrows)
                : make_partition(nrows, (ipart == 1) ? 0 : nprocs-1, ipart);
            AlgPartition cpart = (ipart == 0) ? AlgPartition(ncols)
                : make_partition(ncols, (ipart == 1) ? nprocs-1 : 0, ipart+3);
            auto Ir = IdentityMatrix<Real>(rpart);
            auto Ic = IdentityMatrix<Real>(cpart);

            Real lambda = 2.8;
            int nnz_per_row_max = int(nrows/3);
            auto A = RandomMatrix<Real>(rpart, nrows, ncols, lambda, nnz_per_row_max);

            auto A2 = amrex::SpGEMM(Ir, A, cpart);
            AMREX_ALWAYS_ASSERT(all_true(amrex::almostEqual(A,A2)));

            auto A3 = amrex::SpGEMM(A2, Ic, cpart);
            AMREX_ALWAYS_ASSERT(all_true(amrex::almostEqual(A,A3)));

            auto II = amrex::SpGEMM(Ir, Ir, rpart);
            AMREX_ALWAYS_ASSERT(all_true(amrex::almostEqual(II, Ir)));
        }

        // Permutation Matrix
        {
            std::random_device rd;
            std::uniform_int_distribution<unsigned> dist(0, std::numeric_limits<unsigned>::max());
            unsigned seed = dist(rd);
            ParallelDescriptor::Bcast(&seed, 1);
            std::mt19937 gen(seed);

            int nrows = 240;
            Gpu::PinnedVector<Long> perm(nrows);
            std::iota(perm.begin(), perm.end(), 0);
            std::shuffle(perm.begin(), perm.end(), gen);

            Gpu::DeviceVector<Long> perm_dv(nrows);
            Gpu::copyAsync(Gpu::hostToDevice, perm.begin(), perm.end(),
                           perm_dv.begin());
            auto const* pp = perm_dv.data();

            SpMatrix<float> P(make_partition(nrows, -1, 1), 1);
            P.setVal([=] AMREX_GPU_DEVICE (Long row, Long* col, float* val)
            {
                *col = pp[row];
                *val = 1.0F;
            }, CsrSorted{true});

            auto PT = amrex::transpose(P, P.partition());

            auto PPT = amrex::SpGEMM(P, PT, P.partition());
            auto PTP = amrex::SpGEMM(PT, P, P.partition());
            auto I = IdentityMatrix<float>(P.partition());

            AMREX_ALWAYS_ASSERT(all_true(amrex::almostEqual(PPT,PTP)));
            AMREX_ALWAYS_ASSERT(all_true(amrex::almostEqual(PPT,I)));
        }

        // 1D periodic Laplacian: Lap*Lap == Lap2
        for (int nrows : {6, 128})
        {
            AlgPartition part = (nrows == 6) ? AlgPartition(nrows)
                                             : make_partition(nrows, -1, 2);
            SpMatrix<Real> Lap(part, 3);
            Lap.setVal([=] AMREX_GPU_DEVICE (Long row, Long* col, Real* val)
            {
                if (row == 0) {
                    col[0] = 0;
                    col[1] = 1;
                    col[2] = nrows-1;
                    val[0] = Real(2);
                    val[1] = Real(-1);
                    val[2] = Real(-1);
                } else if (row < nrows-1) {
                    col[0] = row-1;
                    col[1] = row;
                    col[2] = row+1;
                    val[0] = Real(-1);
                    val[1] = Real(2);
                    val[2] = Real(-1);
                } else {
                    col[0] = 0;
                    col[1] = row-1;
                    col[2] = row;
                    val[0] = Real(-1);
                    val[1] = Real(-1);
                    val[2] = Real(2);
                }
            }, CsrSorted{true});

            SpMatrix<Real> Lap2(Lap.partition(), 5);
            Lap2.setVal([=] AMREX_GPU_DEVICE (Long row, Long* col, Real* val)
            {
                if (row == 0) {
                    col[0] = 0;
                    col[1] = 1;
                    col[2] = 2;
                    col[3] = nrows-2;
                    col[4] = nrows-1;
                    val[0] = Real(6);
                    val[1] = Real(-4);
                    val[2] = Real(1);
                    val[3] = Real(1);
                    val[4] = Real(-4);
                } else if (row == 1) {
                    col[0] = 0;
                    col[1] = 1;
                    col[2] = 2;
                    col[3] = 3;
                    col[4] = nrows-1;
                    val[0] = Real(-4);
                    val[1] = Real(6);
                    val[2] = Real(-4);
                    val[3] = Real(1);
                    val[4] = Real(1);
                } else if (row < nrows-2) {
                    col[0] = row-2;
                    col[1] = row-1;
                    col[2] = row;
                    col[3] = row+1;
                    col[4] = row+2;
                    val[0] = Real(1);
                    val[1] = Real(-4);
                    val[2] = Real(6);
                    val[3] = Real(-4);
                    val[4] = Real(1);
                } else if (row == nrows-2) {
                    col[0] = 0;
                    col[1] = row-2;
                    col[2] = row-1;
                    col[3] = row;
                    col[4] = row+1;
                    val[0] = Real(1);
                    val[1] = Real(1);
                    val[2] = Real(-4);
                    val[3] = Real(6);
                    val[4] = Real(-4);
                } else { // row == nrows-1
                    col[0] = 0;
                    col[1] = 1;
                    col[2] = row-2;
                    col[3] = row-1;
                    col[4] = row;
                    val[0] = Real(-4);
                    val[1] = Real(1);
                    val[2] = Real(1);
                    val[3] = Real(-4);
                    val[4] = Real(6);
                }
            }, CsrSorted{true});

            auto LL = amrex::SpGEMM(Lap, Lap, Lap.partition());
            AMREX_ALWAYS_ASSERT(all_true(amrex::almostEqual(LL,Lap2)));
        }

        // (A*B)^T = B^T * A^T, and (A*B)x = A(Bx)
        for (int ipart = 0; ipart < 2; ++ipart)
        {
            Long n1 = 75;
            Long n2 = 100;
            Long n3 = 80;
            AlgPartition pt1 = (ipart == 0) ? AlgPartition(n1) : make_partition(n1, 0, 4);
            AlgPartition pt2 = (ipart == 0) ? AlgPartition(n2) : make_partition(n2, nprocs/2, 5);
            AlgPartition pt3 = (ipart == 0) ? AlgPartition(n3) : make_partition(n3, nprocs-1, 6);
            Real lambda = 3.4;
            int nnz_per_row_max = 9;
            auto A = RandomMatrix<Real>(pt1, n1, n2, lambda, nnz_per_row_max);
            auto B = RandomMatrix<Real>(pt2, n2, n3, lambda, nnz_per_row_max);
            auto AT = amrex::transpose(A, pt2);
            auto BT = amrex::transpose(B, pt3);
            auto AB = amrex::SpGEMM(A, B, pt3);
            auto ABT = amrex::transpose(AB, pt3);
            auto BTAT = amrex::SpGEMM(BT, AT, pt1);
            // Products are summed in different orders.
            AMREX_ALWAYS_ASSERT(all_true(amrex::almostEqual(ABT,BTAT,32)));

            Real tol = std::numeric_limits<Real>::epsilon() * Real(100);
            AMREX_ALWAYS_ASSERT(all_true(check_spmv(A, B, AB, pt3, tol)));
        }

        // Empty matrices
        {
            Long n1 = 40, n2 = 50, n3 = 45;
            AlgPartition pt1 = make_partition(n1, -1, 7);
            AlgPartition pt2 = make_partition(n2, -1, 8);
            AlgPartition pt3 = make_partition(n3, -1, 9);
            SpMatrix<Real> Z1(pt1, 0);
            SpMatrix<Real> Z2(pt2, 0);
            auto A = RandomMatrix<Real>(pt1, n1, n2, 3.0, 6);
            auto B = RandomMatrix<Real>(pt2, n2, n3, 3.0, 6);
            auto ZB = amrex::SpGEMM(Z1, B, pt3);
            auto AZ = amrex::SpGEMM(A, Z2, pt3);
            AMREX_ALWAYS_ASSERT(all_true(ZB.numLocalNonZeros() == 0 &&
                                         AZ.numLocalNonZeros() == 0 &&
                                         ZB.numLocalRows() == pt1.numLocalRows() &&
                                         AZ.numLocalRows() == pt1.numLocalRows()));
            AlgVector<Real> x(pt3), y(pt1);
            x.setVal(Real(1));
            SpMV(y, ZB, x);
            SpMV(y, AZ, x);
            AMREX_ALWAYS_ASSERT(all_true(y.norminf() == Real(0)));
        }
    }
    amrex::Finalize();
}
