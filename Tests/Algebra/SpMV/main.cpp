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
    auto *matCols = mat.columnIndex();

    // simple algebraic system:
    //
    //     s_n + 2.0 * s_{n+1} = (n + 1) + 2 * (n + 2)
    //
    // independent on each processor
    ParallelFor(nrows, [=] AMREX_GPU_DEVICE(Long lrow) {
      auto row = lrow + ib; // global row index

      if (lrow == 0) {
        rhs[0] = Real(1.0);
        phi[0] = Real(1.0);

        matCols[0] = ib;
        matCols[1] = ib + 1;
        matVals[0] = Real(1.0);
        matVals[1] = Real(0.0);
      } else {
        rhs[lrow] = static_cast<Real>(lrow + 2 * (lrow + 1));
        phi[lrow] = static_cast<Real>(lrow + 1);

        matCols[2*lrow    ] = row;
        matCols[2*lrow + 1] = row - 1;
        matVals[2*lrow    ] = static_cast<Real>(2.0);
        matVals[2*lrow + 1] = static_cast<Real>(1.0);
      }
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
