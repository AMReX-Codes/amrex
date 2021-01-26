#include "MyTest.H"
#include <AMReX_Config.H>
#include <AMReX_EB2.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

#if (AMREX_SPACEDIM == 2)
void MyTest::initializePoiseuilleDataFor2D(int ilev) {
  const auto dx = geom[ilev].CellSizeArray();
  for (MFIter mfi(phi[ilev]); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.fabbox();
    Array4<Real> const &fab = phi[ilev].array(mfi);
    Array4<Real> const &fab_gx = grad_x_analytic[ilev].array(mfi);
    Array4<Real> const &fab_gy = grad_y_analytic[ilev].array(mfi);
    Array4<Real> const &fab_eb = grad_eb_analytic[ilev].array(mfi);

    const FabArray<EBCellFlagFab> *flags =
        &(factory[ilev]->getMultiEBCellFlagFab());
    Array4<EBCellFlag const> const &flag = flags->const_array(mfi);

    Array4<Real const> const &ccent = (factory[ilev]->getCentroid()).array(mfi);
    Array4<Real const> const &fcx =
        (factory[ilev]->getFaceCent())[0]->const_array(mfi);
    Array4<Real const> const &fcy =
        (factory[ilev]->getFaceCent())[1]->const_array(mfi);
    Array4<Real const> const &apx =
        (factory[ilev]->getAreaFrac())[0]->const_array(mfi);
    Array4<Real const> const &apy =
        (factory[ilev]->getAreaFrac())[1]->const_array(mfi);
    Array4<Real const> const &norm =
        (factory[ilev]->getBndryNormal()).array(mfi);
    Array4<Real const> const &bcent =
        (factory[ilev]->getBndryCent()).array(mfi);

    const auto &dlo = geom[ilev].Domain().loVect();
    const auto &dhi = geom[ilev].Domain().hiVect();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      Real H = poiseuille_1d_height;
      Real t = (poiseuille_1d_rotation / 180.) * M_PI;

      Real a = std::tan(t);
      Real b = -1.0;
      Real c = poiseuille_1d_pt_on_top_wall[1] -
               std::tan(t) * poiseuille_1d_pt_on_top_wall[0];

      Real rx = (i + 0.5 + ccent(i, j, k, 0)) * dx[0];
      Real ry = (j + 0.5 + ccent(i, j, k, 1)) * dx[1];

      // if not periodic, set the ghost cell values to corr. domain face values
      if (i < dlo[0] and not is_periodic[0]) {
        rx = dlo[0] * dx[0];
        ry = (j + 0.5 + fcx(i, j, k, 0)) * dx[1];
      }
      if (i > dhi[0] and not is_periodic[0]) {
        rx = (dhi[0] + 1) * dx[0];
        ry = (j + 0.5 + fcx(i, j, k, 0)) * dx[1];
      }
      if (j < dlo[1] and not is_periodic[1]) {
        rx = (i + 0.5 + fcy(i, j, k, 0)) * dx[0];
        ry = dlo[1] * dx[1];
      }
      if (j > dhi[1] and not is_periodic[1]) {
        rx = (i + 0.5 + fcy(i, j, k, 0)) * dx[0];
        ry = (dhi[1] + 1) * dx[1];
      }

      auto d = std::fabs(a * rx + b * ry + c) / std::sqrt(a * a + b * b);

//      auto phi_mag = (!flag(i, j, k).isCovered()) ? d * (H - d) : 0.0;
      auto phi_mag = (!flag(i, j, k).isCovered()) ? (H - d) : 0.0;
      fab(i, j, k, 0) = phi_mag * std::cos(t);
      fab(i, j, k, 1) = phi_mag * std::sin(t);

      if (flag(i, j, k).isCovered()) {
        fab_gx(i, j, k, 0) = 0.0;
        fab_gx(i, j, k, 1) = 0.0;
        fab_gy(i, j, k, 0) = 0.0;
        fab_gy(i, j, k, 1) = 0.0;
      } else {

        Real rxl = ((i + 0.5 + ccent(i, j, k, 0))) * dx[0];
        Real ryl = ((j + 0.5 + ccent(i, j, k, 1))) * dx[1];
        Real fac =
          -1.0 ;
        fab_gx(i, j, k, 0) =
                (a * std::cos(t) / std::sqrt(a * a + b * b)) * fac * dx[0];
        fab_gx(i, j, k, 1) =
                 (a * std::sin(t) / std::sqrt(a * a + b * b)) * fac * dx[0];

        rxl = ((i + 0.5 + ccent(i, j, k, 0))) * dx[0];
        ryl = ((j + 0.5 + ccent(i, j, k, 1))) * dx[1];
        fac = -1.0 ;
        fab_gy(i, j, k, 0) =
                (b * std::cos(t) / std::sqrt(a * a + b * b)) * fac * dx[1];
        fab_gy(i, j, k, 1) =
                 (b * std::sin(t) / std::sqrt(a * a + b * b)) * fac * dx[1];
      }

      if (flag(i, j, k).isSingleValued()) {
        Real rxeb = (i + 0.5 + bcent(i, j, k, 0)) * dx[0];
        Real ryeb = (j + 0.5 + bcent(i, j, k, 1)) * dx[1];
        Real fac =
            (H - 2 * (a * rxeb + b * ryeb + c) / (std::sqrt(a * a + b * b)));
        Real dudx = (a * std::cos(t) / std::sqrt(a * a + b * b)) * fac * dx[0];
        Real dvdx = (a * std::sin(t) / std::sqrt(a * a + b * b)) * fac * dx[0];
        Real dudy = (b * std::cos(t) / std::sqrt(a * a + b * b)) * fac * dx[1];
        Real dvdy = (b * std::sin(t) / std::sqrt(a * a + b * b)) * fac * dx[1];

        fab_eb(i, j, k, 0) = dudx * norm(i, j, k, 0) + dudy * norm(i, j, k, 1);
        fab_eb(i, j, k, 1) = dvdx * norm(i, j, k, 0) + dvdy * norm(i, j, k, 1);
      }
    });
  }
}

#else
void MyTest::initializePoiseuilleDataFor2D(int ilev) {
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0,
                                   "Calling 2D function for 3D case. Error!!");
}
#endif
