#include "MyTest.H"
#include <AMReX_EB2.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

#if (AMREX_SPACEDIM == 2)
void MyTest::initializeTriangleWaveDataFor2D(int ilev) {
  const auto dx = geom[ilev].CellSizeArray();
  for (MFIter mfi(phi[ilev]); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.fabbox();
    Array4<Real> const &fab = phi[ilev].array(mfi);
    Array4<Real> const &fabphieb = phieb[ilev].array(mfi);
    Array4<Real> const &fab_ghost_resolved =
        phi_ghost_resolved[ilev].array(mfi);
    Array4<Real> const &fab_gx = grad_x_analytic[ilev].array(mfi);
    Array4<Real> const &fab_gy = grad_y_analytic[ilev].array(mfi);
    Array4<Real> const &fab_eb = grad_eb_analytic[ilev].array(mfi);
    Array4<Real> const &fab_lap = lap_analytic[ilev].array(mfi);

    const FabArray<EBCellFlagFab> *flags =
        &(factory[ilev]->getMultiEBCellFlagFab());
    Array4<EBCellFlag const> const &flag = flags->const_array(mfi);

    const auto &dlo = geom[ilev].Domain().loVect();
    const auto &dhi = geom[ilev].Domain().hiVect();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // y components are all zero
      fab(i, j, k, 1) = 0.0;
      fab_ghost_resolved(i, j, k, 1) = 0.0;
      fabphieb(i, j, k, 1) = 0.0;

      if (flag(i, j, k).isCovered()) {
        fab(i, j, k, 0) = 0.0;
        fab_ghost_resolved(i, j, k, 0) = 0.0;
        fabphieb(i, j, k, 0) = 0.0;
      } else {
        if (std::abs(i) % 2 == 0) {
          fab(i, j, k, 0) = 1.0;
          fabphieb(i, j, k, 0) = 1.0;
          fab_ghost_resolved(i, j, k, 0) = 1.0;
        } else {
          fab(i, j, k, 0) = 2.0;
          fabphieb(i, j, k, 0) = 2.0;
          fab_ghost_resolved(i, j, k, 0) = 2.0;
        }
      }

      // if outside domain and not periodic, phi must be the domain face value
      if (!is_periodic[0] && (i < dlo[0] || i > dhi[0])) {
        fab(i, j, k, 0) = 1.5;
      }

      fab_gx(i, j, k, 1) = 0.0;
      fab_gy(i, j, k, 0) = 0.0;
      fab_gy(i, j, k, 1) = 0.0;
      fab_lap(i, j, k, 1) = 0.0;

      if (flag(i, j, k).isCovered()) {
        fab_gx(i, j, k, 0) = 0.0;
        fab_lap(i, j, k, 0) = 0.0;
      } else {
        if (std::abs(i) % 2 == 0) {
          fab_gx(i, j, k, 0) =
              (-1.0 / dx[0]) * dx[0]; // needs to be multiplied by dx for
                                      // comparison with grad functions
          fab_lap(i, j, k, 0) = -2.0 / (dx[0] * dx[0]);
        } else {
          fab_gx(i, j, k, 0) =
              (1.0 / dx[0]) * dx[0]; // needs to be multiplied by dx for
                                     // comparison with grad functions
          fab_lap(i, j, k, 0) = 2.0 / (dx[0] * dx[0]);
        }

        if (flag(i, j, k).isSingleValued()) {
          fab_eb(i, j, k, 0) = 0.0;
          fab_eb(i, j, k, 1) = 0.0;
        }
      }
    });
  }
}
#endif
