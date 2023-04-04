#include "MyTest.H"
#include <AMReX_EB2.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

#if (AMREX_SPACEDIM == 3)
void MyTest::initializePoiseuilleDataFor3D(int ilev) {
  const auto dx = geom[ilev].CellSizeArray();
  for (MFIter mfi(phi[ilev]); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.fabbox();
    Array4<Real> const &fab = phi[ilev].array(mfi);
    Array4<Real> const &fab_ghost_resolved = phi_ghost_resolved[ilev].array(mfi);
    Array4<Real> const &fab_gx = grad_x_analytic[ilev].array(mfi);
    Array4<Real> const &fab_gy = grad_y_analytic[ilev].array(mfi);
    Array4<Real> const &fab_gz = grad_z_analytic[ilev].array(mfi);
    Array4<Real> const &fab_eb = grad_eb_analytic[ilev].array(mfi);
    Array4<Real> const &fab_lap = lap_analytic[ilev].array(mfi);

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

    Array4<Real const> const &fcz =
        (factory[ilev]->getFaceCent())[2]->const_array(mfi);
    Array4<Real const> const &apz =
        (factory[ilev]->getAreaFrac())[2]->const_array(mfi);

    if (poiseuille_askew) { // 3D askew
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j,
                                                  int k) noexcept {
        Real H = poiseuille_height;
        int nfdir = poiseuille_no_flow_dir;
        constexpr Real pi = 3.1415926535897932;
        Real alpha = (poiseuille_askew_rotation[0] / 180.) * pi;
        Real gamma = (poiseuille_askew_rotation[1] / 180.) * pi;

        Real a = std::sin(gamma);
        Real b = -std::cos(alpha) * std::cos(gamma);
        Real c = std::sin(alpha);
        Real d = -a * poiseuille_pt_on_top_wall[0] -
                 b * poiseuille_pt_on_top_wall[1] -
                 c * poiseuille_pt_on_top_wall[2];

        Real rx = (i + 0.5 + ccent(i, j, k, 0)) * dx[0];
        Real ry = (j + 0.5 + ccent(i, j, k, 1)) * dx[1];
        Real rz = (k + 0.5 + ccent(i, j, k, 2)) * dx[2];
        Real rx_gr = rx;
        Real ry_gr = ry;
        Real rz_gr = rz;

        // if not periodic, set the ghost cell values to corr. domain face
        // values
        if (i < dlo[0] && !is_periodic[0]) {
          rx = dlo[0] * dx[0];
          ry = (j + 0.5 + fcx(i, j, k, 0)) * dx[1];
          rz = (k + 0.5 + fcx(i, j, k, 1)) * dx[2];
        }
        if (i > dhi[0] && !is_periodic[0]) {
          rx = (dhi[0] + 1) * dx[0];
          ry = (j + 0.5 + fcx(i, j, k, 0)) * dx[1];
          rz = (k + 0.5 + fcx(i, j, k, 1)) * dx[2];
        }
        if (j < dlo[1] && !is_periodic[1]) {
          rx = (i + 0.5 + fcy(i, j, k, 0)) * dx[0];
          ry = dlo[1] * dx[1];
          rz = (k + 0.5 + fcy(i, j, k, 1)) * dx[2];
        }
        if (j > dhi[1] && !is_periodic[1]) {
          rx = (i + 0.5 + fcy(i, j, k, 0)) * dx[0];
          ry = (dhi[1] + 1) * dx[1];
          rz = (k + 0.5 + fcy(i, j, k, 1)) * dx[2];
        }
        if (k < dlo[2] && !is_periodic[2]) {
          rx = (i + 0.5 + fcz(i, j, k, 0)) * dx[0];
          ry = (j + 0.5 + fcz(i, j, k, 1)) * dx[1];
          rz = dlo[2] * dx[2];
        }
        if (k > dhi[2] && !is_periodic[2]) {
          rx = (i + 0.5 + fcz(i, j, k, 0)) * dx[0];
          ry = (j + 0.5 + fcz(i, j, k, 1)) * dx[1];
          rz = (dhi[2] + 1) * dx[2];
        }

        auto dist = std::fabs(a * rx + b * ry + c * rz + d) /
                    std::sqrt(a * a + b * b + c * c);

        auto phi_mag = (!flag(i, j, k).isCovered()) ? dist * (H - dist) : 0.0;

        auto dist_gr = std::fabs(a * rx_gr + b * ry_gr + c * rz_gr + d) /
                    std::sqrt(a * a + b * b + c * c);
        auto phi_mag_gr = (!flag(i, j, k).isCovered()) ? dist_gr * (H - dist_gr) : 0.0;

        Vector<Real> flow_norm(3, 0.0);

        if (nfdir == 2) {
          Real flow_norm_mag = std::sqrt(std::cos(alpha) * std::cos(alpha) *
                                             std::cos(gamma) * std::cos(gamma) +
                                         std::sin(gamma) * std::sin(gamma));
          flow_norm[0] = std::cos(alpha) * std::cos(gamma) / flow_norm_mag;
          flow_norm[1] = std::sin(gamma) / flow_norm_mag;
        } else if (nfdir == 1) {
          Real flow_norm_mag = std::sqrt(std::sin(alpha) * std::sin(alpha) +
                                         std::sin(gamma) * std::sin(gamma));
          flow_norm[0] = -std::sin(alpha) / flow_norm_mag;
          flow_norm[2] = std::sin(gamma) / flow_norm_mag;

        } else if (nfdir == 0) {
          Real flow_norm_mag = std::sqrt(std::cos(alpha) * std::cos(alpha) *
                                             std::cos(gamma) * std::cos(gamma) +
                                         std::sin(alpha) * std::sin(alpha));
          flow_norm[2] = std::cos(alpha) * std::cos(gamma) / flow_norm_mag;
          flow_norm[1] = std::sin(alpha) / flow_norm_mag;
        } else {
          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(1 == 1, "Invalid flow direction");
        }

        fab(i, j, k, 0) = phi_mag * flow_norm[0];
        fab(i, j, k, 1) = phi_mag * flow_norm[1];
        fab(i, j, k, 2) = phi_mag * flow_norm[2];

        fab_ghost_resolved(i, j, k, 0) = phi_mag_gr * flow_norm[0];
        fab_ghost_resolved(i, j, k, 1) = phi_mag_gr * flow_norm[1];
        fab_ghost_resolved(i, j, k, 2) = phi_mag_gr * flow_norm[2];

        if (flag(i, j, k).isCovered()) {
          fab_gx(i, j, k, 0) = 0.0;
          fab_gx(i, j, k, 1) = 0.0;
          fab_gx(i, j, k, 2) = 0.0;
          fab_gy(i, j, k, 0) = 0.0;
          fab_gy(i, j, k, 1) = 0.0;
          fab_gy(i, j, k, 2) = 0.0;
          fab_gz(i, j, k, 0) = 0.0;
          fab_gz(i, j, k, 1) = 0.0;
          fab_gz(i, j, k, 2) = 0.0;

          fab_lap(i, j, k, 0) = 0.0;
          fab_lap(i, j, k, 1) = 0.0;
          fab_lap(i, j, k, 2) = 0.0;
        } else {
          fab_lap(i, j, k, 0) =
              2.0 * a * a * flow_norm[0] / (a * a + b * b + c * c) +
              2.0 * b * b * flow_norm[0] / (a * a + b * b + c * c) +
              2.0 * c * c * flow_norm[0] / (a * a + b * b + c * c);
          fab_lap(i, j, k, 1) =
              2.0 * a * a * flow_norm[1] / (a * a + b * b + c * c) +
              2.0 * b * b * flow_norm[1] / (a * a + b * b + c * c) +
              2.0 * c * c * flow_norm[1] / (a * a + b * b + c * c);
          fab_lap(i, j, k, 2) =
              2.0 * a * a * flow_norm[2] / (a * a + b * b + c * c) +
              2.0 * b * b * flow_norm[2] / (a * a + b * b + c * c) +
              2.0 * c * c * flow_norm[2] / (a * a + b * b + c * c);

          Real rxl = i * dx[0];
          Real ryl = (j + 0.5 + fcx(i, j, k, 0)) * dx[1];
          Real rzl = (k + 0.5 + fcx(i, j, k, 1)) * dx[2];
          Real fac = (H - 2 * (a * rxl + b * ryl + c * rzl + d) /
                              (std::sqrt(a * a + b * b + c * c)));
          fab_gx(i, j, k, 0) =
              (apx(i, j, k) == 0.0)
                  ? 0.0
                  : (a * flow_norm[0] / std::sqrt(a * a + b * b + c * c)) *
                        fac * dx[0];
          fab_gx(i, j, k, 1) =
              (apx(i, j, k) == 0.0)
                  ? 0.0
                  : (a * flow_norm[1] / std::sqrt(a * a + b * b + c * c)) *
                        fac * dx[0];
          fab_gx(i, j, k, 2) =
              (apx(i, j, k) == 0.0)
                  ? 0.0
                  : (a * flow_norm[2] / std::sqrt(a * a + b * b + c * c)) *
                        fac * dx[0];

          rxl = (i + 0.5 + fcy(i, j, k, 0)) * dx[0];
          ryl = j * dx[1];
          rzl = (k + 0.5 + fcy(i, j, k, 1)) * dx[2];
          fac = (H - 2 * (a * rxl + b * ryl + c * rzl + d) /
                         (std::sqrt(a * a + b * b + c * c)));
          fab_gy(i, j, k, 0) =
              (apy(i, j, k) == 0.0)
                  ? 0.0
                  : (b * flow_norm[0] / std::sqrt(a * a + b * b + c * c)) *
                        fac * dx[1];
          fab_gy(i, j, k, 1) =
              (apy(i, j, k) == 0.0)
                  ? 0.0
                  : (b * flow_norm[1] / std::sqrt(a * a + b * b + c * c)) *
                        fac * dx[1];
          fab_gy(i, j, k, 2) =
              (apy(i, j, k) == 0.0)
                  ? 0.0
                  : (b * flow_norm[2] / std::sqrt(a * a + b * b + c * c)) *
                        fac * dx[1];

          rxl = (i + 0.5 + fcz(i, j, k, 0)) * dx[0];
          ryl = (j + 0.5 + fcz(i, j, k, 1)) * dx[1];
          rzl = k * dx[2];
          fac = (H - 2 * (a * rxl + b * ryl + c * rzl + d) /
                         (std::sqrt(a * a + b * b + c * c)));
          fab_gz(i, j, k, 0) =
              (apz(i, j, k) == 0.0)
                  ? 0.0
                  : (c * flow_norm[0] / std::sqrt(a * a + b * b + c * c)) *
                        fac * dx[2];
          fab_gz(i, j, k, 1) =
              (apz(i, j, k) == 0.0)
                  ? 0.0
                  : (c * flow_norm[1] / std::sqrt(a * a + b * b + c * c)) *
                        fac * dx[2];
          fab_gz(i, j, k, 2) =
              (apz(i, j, k) == 0.0)
                  ? 0.0
                  : (c * flow_norm[2] / std::sqrt(a * a + b * b + c * c)) *
                        fac * dx[2];
        }

        if (flag(i, j, k).isSingleValued()) {
          Real rxeb = (i + 0.5 + bcent(i, j, k, 0)) * dx[0];
          Real ryeb = (j + 0.5 + bcent(i, j, k, 1)) * dx[1];
          Real rzeb = (k + 0.5 + bcent(i, j, k, 2)) * dx[2];
          Real fac = (H - 2 * (a * rxeb + b * ryeb + c * rzeb + d) /
                              (std::sqrt(a * a + b * b + c * c)));
          Real dudx = (a * flow_norm[0] / std::sqrt(a * a + b * b + c * c)) *
                      fac * dx[0];
          Real dvdx = (a * flow_norm[1] / std::sqrt(a * a + b * b + c * c)) *
                      fac * dx[0];
          Real dwdx = (a * flow_norm[2] / std::sqrt(a * a + b * b + c * c)) *
                      fac * dx[0];
          Real dudy = (b * flow_norm[0] / std::sqrt(a * a + b * b + c * c)) *
                      fac * dx[1];
          Real dvdy = (b * flow_norm[1] / std::sqrt(a * a + b * b + c * c)) *
                      fac * dx[1];
          Real dwdy = (b * flow_norm[2] / std::sqrt(a * a + b * b + c * c)) *
                      fac * dx[1];
          Real dudz = (c * flow_norm[0] / std::sqrt(a * a + b * b + c * c)) *
                      fac * dx[2];
          Real dvdz = (c * flow_norm[1] / std::sqrt(a * a + b * b + c * c)) *
                      fac * dx[2];
          Real dwdz = (c * flow_norm[2] / std::sqrt(a * a + b * b + c * c)) *
                      fac * dx[2];

          fab_eb(i, j, k, 0) = dudx * norm(i, j, k, 0) +
                               dudy * norm(i, j, k, 1) +
                               dudz * norm(i, j, k, 2);
          fab_eb(i, j, k, 1) = dvdx * norm(i, j, k, 0) +
                               dvdy * norm(i, j, k, 1) +
                               dvdz * norm(i, j, k, 2);
          fab_eb(i, j, k, 2) = dwdx * norm(i, j, k, 0) +
                               dwdy * norm(i, j, k, 1) +
                               dwdz * norm(i, j, k, 2);
        }

      });
    } else { // 3D grid-aligned
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j,
                                                  int k) noexcept {
        Real H = poiseuille_height;
        Real bot = poiseuille_bottom;
        int dir = poiseuille_height_dir;
        int fdir = poiseuille_flow_dir;

        fab(i, j, k, 0) = 0.0;
        fab(i, j, k, 1) = 0.0;
        fab(i, j, k, 2) = 0.0;
        fab_gx(i, j, k, 0) = 0.0;
        fab_gx(i, j, k, 1) = 0.0;
        fab_gx(i, j, k, 2) = 0.0;
        fab_gy(i, j, k, 0) = 0.0;
        fab_gy(i, j, k, 1) = 0.0;
        fab_gy(i, j, k, 2) = 0.0;
        fab_gz(i, j, k, 0) = 0.0;
        fab_gz(i, j, k, 1) = 0.0;
        fab_gz(i, j, k, 2) = 0.0;
        fab_lap(i, j, k, 0) = 0.0;
        fab_lap(i, j, k, 1) = 0.0;
        fab_lap(i, j, k, 2) = 0.0;

        Real d = 0.0;
        Real d_gr = 0.0;
        if (dir == 0) {
          Real rx = (i + 0.5 + ccent(i, j, k, 0)) * dx[0];
          Real rx_gr = rx;

          if (i < dlo[0] && !is_periodic[0]) {
            rx = dlo[0] * dx[0];
          }
          if (i > dhi[0] && !is_periodic[0]) {
            rx = (dhi[0] + 1) * dx[0];
          }

          d = rx - bot;
          d_gr = rx_gr - bot;

          fab(i, j, k, fdir) = (!flag(i, j, k).isCovered()) ? d * (H - d) : 0.0;
          fab_lap(i, j, k, fdir) = (!flag(i, j, k).isCovered()) ? 2.0 : 0.0;

          fab_ghost_resolved(i, j, k, fdir) = (!flag(i, j, k).isCovered()) ? d_gr * (H - d_gr) : 0.0;

          Real rxl = i * dx[0];
          d = rxl - bot;
          fab_gx(i, j, k, fdir) =
              (apx(i, j, k) == 0.0) ? 0.0 : (H - 2 * d) * dx[0];

          if (flag(i, j, k).isSingleValued()) {
            fab_eb(i, j, k, 0) = 0.0;
            fab_eb(i, j, k, 1) = 0.0;
            fab_eb(i, j, k, 2) = 0.0;

            Real rxeb = (i + 0.5 + bcent(i, j, k, 0)) * dx[0];
            d = rxeb - bot;
            fab_eb(i, j, k, fdir) = (H - 2 * d) * dx[0] * norm(i, j, k, 0);
          }
        } else if (dir == 1) {
          Real ry = (j + 0.5 + ccent(i, j, k, 1)) * dx[1];
          Real ry_gr = ry;

          if (j < dlo[1] && !is_periodic[1]) {
            ry = dlo[1] * dx[1];
          }
          if (j > dhi[1] && !is_periodic[1]) {
            ry = (dhi[1] + 1) * dx[1];
          }

          d = ry - bot;
          d_gr = ry_gr - bot;

          fab(i, j, k, fdir) = (!flag(i, j, k).isCovered()) ? d * (H - d) : 0.0;
          fab_lap(i, j, k, fdir) = (!flag(i, j, k).isCovered()) ? 2.0 : 0.0;

          fab_ghost_resolved(i, j, k, fdir) = (!flag(i, j, k).isCovered()) ? d_gr * (H - d_gr) : 0.0;

          Real ryl = j * dx[1];
          d = ryl - bot;
          fab_gy(i, j, k, fdir) =
              (apy(i, j, k) == 0.0) ? 0.0 : (H - 2 * d) * dx[1];

          if (flag(i, j, k).isSingleValued()) {
            fab_eb(i, j, k, 0) = 0.0;
            fab_eb(i, j, k, 1) = 0.0;
            fab_eb(i, j, k, 2) = 0.0;

            Real ryeb = (j + 0.5 + bcent(i, j, k, 1)) * dx[1];
            d = ryeb - bot;
            fab_eb(i, j, k, fdir) = (H - 2 * d) * dx[1] * norm(i, j, k, 1);
          }
        } else if (dir == 2) {
          Real rz = (k + 0.5 + ccent(i, j, k, 2)) * dx[2];
          Real rz_gr = rz;

          if (k < dlo[2] && !is_periodic[2]) {
            rz = dlo[2] * dx[2];
          }
          if (k > dhi[2] && !is_periodic[2]) {
            rz = (dhi[2] + 1) * dx[2];
          }

          d = rz - bot;
          d_gr = rz_gr - bot;

          fab(i, j, k, fdir) = (!flag(i, j, k).isCovered()) ? d * (H - d) : 0.0;
          fab_lap(i, j, k, fdir) = (!flag(i, j, k).isCovered()) ? 2.0 : 0.0;

          fab_ghost_resolved(i, j, k, fdir) = (!flag(i, j, k).isCovered()) ? d_gr * (H - d_gr) : 0.0;

          Real rzl = k * dx[2];
          d = rzl - bot;
          fab_gz(i, j, k, fdir) =
              (apz(i, j, k) == 0.0) ? 0.0 : (H - 2 * d) * dx[2];

          if (flag(i, j, k).isSingleValued()) {
            fab_eb(i, j, k, 0) = 0.0;
            fab_eb(i, j, k, 1) = 0.0;
            fab_eb(i, j, k, 2) = 0.0;

            Real rzeb = (k + 0.5 + bcent(i, j, k, 2)) * dx[2];
            d = rzeb - bot;
            fab_eb(i, j, k, fdir) = (H - 2 * d) * dx[2] * norm(i, j, k, 2);
          }
        } else {
          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0, "Invalid height direction");
        }
      });
    }
  }
}

#else
void MyTest::initializePoiseuilleDataFor3D(int ilev) {
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(0,
                                   "Calling 3D function for 2D case. Error!!");
}
#endif
