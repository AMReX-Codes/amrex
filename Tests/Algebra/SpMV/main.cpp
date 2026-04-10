#include <AMReX_GMRES_MV.H>
#include <AMReX_SpMV.H>

#include <AMReX.H>

#include <cmath>

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
                   << multiplicationError << ", solve = " << solveError << "\n\n";

    AMREX_ALWAYS_ASSERT(multiplicationError < eps && solveError < eps);
  }

  // restriction
  {
      amrex::Print() << "Restriction Test:\n";

      using T = amrex::Real;

      int ncells_f = 32;
      BoxND<2> fbox(IntVectND<2>(0), IntVectND<2>(ncells_f-1));
      AlgVector<T> fvec(fbox.numPts());
      auto nrows_f = fvec.numLocalRows();
      auto row_begin_f = fvec.globalBegin();
      auto* pf = fvec.data();

      BoxND<2> cbox = amrex::coarsen(fbox, 2);
      AlgVector<T> cvec(cbox.numPts());
      auto nrows_c = cvec.numLocalRows();
      auto row_begin_c = cvec.globalBegin();
      auto* pc = cvec.data();

      ParallelFor(nrows_f, [=] AMREX_GPU_DEVICE (Long lrow)
      {
          auto row = lrow + row_begin_f;
          auto civ = amrex::coarsen(fbox.atOffset(row),2);
          auto x = civ[0] + 0.5;
          auto y = civ[1] + 0.5;
          pf[lrow] = T(std::hypot(x, y));
      });

      SpMatrix<T> rmat(cvec.partition(), 4);

      rmat.setVal([=] AMREX_GPU_DEVICE (Long row, Long* col, T* val) {
          auto civ = cbox.atOffset(row);
          int count = 0;
          for (int ry = 0; ry < 2; ++ry) {
          for (int rx = 0; rx < 2; ++rx) {
              auto fiv = civ*2 + IntVectND<2>(rx,ry);
              col[count] = fbox.index(fiv);
              val[count] = T(0.25*(1.0 + 1.e-2*(2*rx-1) + 1.e-1*(2*ry-1)));
              ++count;
          }}
      }, CsrSorted{false});

      amrex::SpMV(cvec, rmat, fvec);

      ParallelFor(nrows_c, [=] AMREX_GPU_DEVICE (Long lrow)
      {
          auto row = lrow + row_begin_c;
          auto iv = cbox.atOffset(row);
          auto x = iv[0] + 0.5;
          auto y = iv[1] + 0.5;
          auto expected = T(std::hypot(x, y));
          pc[lrow] -= expected;
      });

      auto error = cvec.norminf();
      amrex::Print() << " Max norm error of restriction: " << error << "\n\n";
      auto eps = (sizeof(T) == 4) ? T(1.e-5) : T(1.e-14);
      AMREX_ALWAYS_ASSERT(error < eps);
  }

  // interpolation
  {
      using T = std::conditional_t<std::is_same_v<Real,double>,float,double>;

      amrex::Print() << "Interpolation Test (" << amrex::demangle(typeid(T).name()) << ")\n";

      int ncells_c = 15;
      BoxND<2> cbox(IntVectND<2>(0), IntVectND<2>(ncells_c), IntVectND<2>(1));
      AlgVector<T> cvec(cbox.numPts());
      auto nrows_c = cvec.numLocalRows();
      auto row_begin_c = cvec.globalBegin();
      auto* pc = cvec.data();

      BoxND<2> fbox = amrex::refine(cbox, 2);
      AlgVector<T> fvec(fbox.numPts());
      auto nrows_f = fvec.numLocalRows();
      auto row_begin_f = fvec.globalBegin();
      auto* pf = fvec.data();

      double dxc = 0.12;
      double dyc = 0.19;

      ParallelFor(nrows_c, [=] AMREX_GPU_DEVICE (Long lrow)
      {
          auto row = lrow + row_begin_c;
          auto civ = cbox.atOffset(row);
          auto x = civ[0] * dxc;
          auto y = civ[1] * dyc;
          pc[lrow] = T(x+y);
      });

      Long nnz_max = nrows_f * 4;
      Long actual_nnz;
      Gpu::DeviceVector<T> mat_dv(nnz_max);
      Gpu::DeviceVector<Long> col_dv(nnz_max,-1);
      Gpu::DeviceVector<Long> row_dv(nrows_f+1);
      {
          auto* mat = mat_dv.data();
          auto* col = col_dv.data();
          ParallelFor(nrows_f, [=] AMREX_GPU_DEVICE (Long lrow)
          {
              auto row = lrow + row_begin_f;
              auto fiv = fbox.atOffset(row);
              auto civ = amrex::coarsen(fiv,2);
              bool x_is_even = civ[0]*2 == fiv[0];
              bool y_is_even = civ[1]*2 == fiv[1];
              if (x_is_even && y_is_even) {
                  mat[4*lrow] = T(1);
                  col[4*lrow] = cbox.index(civ);
              } else if (x_is_even) {
                  mat[4*lrow  ] = T(0.5);
                  mat[4*lrow+1] = T(0.5);
                  col[4*lrow  ] = cbox.index(civ);
                  col[4*lrow+1] = cbox.index(civ + IntVectND<2>{0,1});
              } else if (y_is_even) {
                  mat[4*lrow  ] = T(0.5);
                  mat[4*lrow+1] = T(0.5);
                  col[4*lrow  ] = cbox.index(civ);
                  col[4*lrow+1] = cbox.index(civ + IntVectND<2>{1,0});
              } else {
                  mat[4*lrow  ] = T(0.25);
                  mat[4*lrow+1] = T(0.25);
                  mat[4*lrow+2] = T(0.25);
                  mat[4*lrow+3] = T(0.25);
                  col[4*lrow  ] = cbox.index(civ);
                  col[4*lrow+1] = cbox.index(civ + IntVectND<2>{0,1});
                  col[4*lrow+2] = cbox.index(civ + IntVectND<2>{1,0});
                  col[4*lrow+3] = cbox.index(civ + IntVectND<2>{1,1});
              }
          });

          // For this small problem, let's use int.
          Gpu::DeviceVector<int> psum_dv(nnz_max);
          auto* ps = psum_dv.data();
          actual_nnz = Scan::PrefixSum<int>
              (int(nnz_max),
               [=] AMREX_GPU_DEVICE (int i) -> int { return col[i] >= 0; },
               [=] AMREX_GPU_DEVICE (int i, int x) { ps[i] = x; },
               Scan::Type::exclusive, Scan::retSum);

          Gpu::DeviceVector<T> mat_tmp(actual_nnz);
          Gpu::DeviceVector<Long> col_tmp(actual_nnz);
          auto* mat2 = mat_tmp.data();
          auto* col2 = col_tmp.data();
          ParallelFor(nnz_max, [=] AMREX_GPU_DEVICE (Long i)
          {
              if (col[i] >= 0) {
                  mat2[ps[i]] = mat[i];
                  col2[ps[i]] = col[i];
              }
          });

          auto* row_offset = row_dv.data();
          ParallelFor(nrows_f, [=] AMREX_GPU_DEVICE (Long i)
          {
              row_offset[i] = ps[i*4];
              if (i == nrows_f - 1) {
                  row_offset[nrows_f] = actual_nnz;
              }
          });

          std::swap(mat_dv, mat_tmp);
          std::swap(col_dv, col_tmp);
          Gpu::streamSynchronize();
      }

      SpMatrix<T> pmat{};
      pmat.define(fvec.partition(), mat_dv.data(), col_dv.data(), actual_nnz,
                  row_dv.data(), CsrSorted{false}, CsrValid{true});

      amrex::SpMV(fvec, pmat, cvec);

      ParallelFor(nrows_f, [=] AMREX_GPU_DEVICE (Long lrow)
      {
          auto row = lrow + row_begin_f;
          auto fiv = fbox.atOffset(row);
          auto x = fiv[0] * dxc * 0.5;
          auto y = fiv[1] * dyc * 0.5;
          auto expected = T(x+y);
          pf[lrow] -= expected;
      });

      auto error = fvec.norminf();
      amrex::Print() << " Max norm error of interpolation: " << error << "\n\n";
      auto eps = (sizeof(T) == 4) ? T(1.e-5) : T(1.e-14);
      AMREX_ALWAYS_ASSERT(error < eps);
  }

  amrex::Finalize();
}
