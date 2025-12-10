#include <AMReX_GMRES_MV.H>
#include <AMReX_SpMV.H>

#include <AMReX.H>

using namespace amrex;

int main(int argc, char *argv[]) {
  amrex::Initialize(argc, argv);
  {
    AlgVector<Real> xvec(100);
    AlgVector<Real> bvec(xvec.partition());
    AlgVector<Real> exact(xvec.partition());

    int num_non_zeros = 2;
    SpMatrix<Real> mat(xvec.partition(), num_non_zeros);

    auto *rhs = bvec.data();
    auto *phi = exact.data();
    auto nrows = bvec.numLocalRows();
    auto ib = bvec.globalBegin();

    auto *matVals = mat.data();
    auto *matRowOffsets = mat.rowOffset();
    auto *matCols = mat.columnIndex();

    // simple algebraic system:
    //
    //     s_n + 2.0 * s_{n+1} = (n + 1) + 2 * (n + 2)
    //
    // independent on each processor
    ParallelFor(nrows - 1, [=] AMREX_GPU_DEVICE(Long lrow) {
      auto row = lrow + ib; // global row index

      rhs[lrow + 1] = static_cast<Real>((lrow + 1.0) + 2.0 * (lrow + 2.0));
      phi[lrow + 1] = static_cast<Real>(lrow + 2);

      matRowOffsets[lrow + 1] = static_cast<Long>(2 * (lrow + 1));
      matCols[2 * (lrow + 1)] = row + 1;
      matCols[2 * (lrow + 1) + 1] = row;
      matVals[2 * (lrow + 1)] = static_cast<Real>(2.0);
      matVals[2 * (lrow + 1) + 1] = static_cast<Real>(1.0);
    });

    // boundary condition: s_0 = 1.0
    ParallelFor(1, [=] AMREX_GPU_DEVICE(Long) {
      // global row index

      rhs[0] = 1.0;
      phi[0] = 1.0;

      matRowOffsets[0] = 0;
      matCols[0] = ib;
      matCols[1] = ib + 1;
      matVals[0] = 1.0;
      matVals[1] = 0.0;
    });

    auto eps = (sizeof(Real) == 4) ? Real(1.e-5) : Real(1.e-12);
    amrex::SpMV(xvec, mat, exact);

    // Check the multiplication
    amrex::Axpy(xvec, Real(-1.0), bvec);

    Real multiplicationError = xvec.norminf();

    xvec.setVal(1.0);

    GMRES_MV<Real> gmres(&mat);
    gmres.setPrecond(JacobiSmoother<Real>(&mat));
    gmres.setVerbose(2);

    gmres.solve(xvec, bvec, Real(0.0), eps);

    // Check the solution
    amrex::Axpy(xvec, Real(-1.0), exact);

    auto solveError = xvec.norminf();
    amrex::Print() << " Max norm error: multiplication = "
                   << multiplicationError << ", solve = " << solveError << "\n";

    AMREX_ALWAYS_ASSERT(multiplicationError < eps && solveError < eps);
  }
  amrex::Finalize();
}
