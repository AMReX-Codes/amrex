#ifndef MY_TEST_K_H_
#define MY_TEST_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_EBCellFlag.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mytest_set_phi_reg (int i, int j, int k, amrex::Array4<amrex::Real> const& phie,
                         amrex::Array4<amrex::Real> const& rhs,
                         AMREX_D_DECL(amrex::Array4<amrex::Real> const& bx,
                                      amrex::Array4<amrex::Real> const& by,
                                      amrex::Array4<amrex::Real> const& bz),
                         amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                         int prob_type, amrex::Box const& vbx)
{
    constexpr amrex::Real pi = 3.1415926535897932;
    amrex::IntVect iv(AMREX_D_DECL(i,j,k));

    if (vbx.contains(iv)) {
        amrex::Real x = (i+0.5)*dx[0] - 0.5;
        amrex::Real y = (j+0.5)*dx[1] - 0.5;
        amrex::Real theta = std::atan2(x,y) + 0.5*pi;
        amrex::Real r2 = x*x+y*y;
        phie(i,j,k) = r2*r2 * std::cos(3.*theta);
        if (prob_type == 1) {
            rhs(i,j,k) = -7.0 * r2 * std::cos(3.*theta);
        } else {
            rhs(i,j,k) = -(7.0 * r2 - 15. * r2*r2) * std::cos(3.*theta);
        }
    }

    if (prob_type == 2) {
        if (amrex::surroundingNodes(vbx,0).contains(iv)) {
            amrex::Real x = (i    )*dx[0] - 0.5;
            amrex::Real y = (j+0.5)*dx[1] - 0.5;
            amrex::Real theta = std::atan2(x,y) + 0.5*pi;
            amrex::Real r2 = x*x+y*y;
            bx(i,j,k) = 1.0 - r2;
        }
        if (amrex::surroundingNodes(vbx,1).contains(iv)) {
            amrex::Real x = (i+0.5)*dx[0] - 0.5;
            amrex::Real y = (j    )*dx[1] - 0.5;
            amrex::Real theta = std::atan2(x,y) + 0.5*pi;
            amrex::Real r2 = x*x+y*y;
            by(i,j,k) = 1.0 - r2;
        }
#if (AMREX_SPACEDIM == 3)
        if (amrex::surroundingNodes(vbx,2).contains(iv)) {
            amrex::Real x = (i+0.5)*dx[0] - 0.5;
            amrex::Real y = (j+0.5)*dx[1] - 0.5;
            amrex::Real theta = std::atan2(x,y) + 0.5*pi;
            amrex::Real r2 = x*x+y*y;
            bz(i,j,k) = 1.0 - r2;
        }
#endif
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mytest_set_phi_eb (int i, int j, int k, amrex::Array4<amrex::Real> const& phie,
                        amrex::Array4<amrex::Real> const& phib,
                        amrex::Array4<amrex::Real> const& rhs,
                        AMREX_D_DECL(amrex::Array4<amrex::Real> const& bx,
                                     amrex::Array4<amrex::Real> const& by,
                                     amrex::Array4<amrex::Real> const& bz),
                        amrex::Array4<amrex::Real> const& bb,
                        amrex::Array4<amrex::EBCellFlag const> const& flag,
                        amrex::Array4<amrex::Real const> const& cent,
                        amrex::Array4<amrex::Real const> const& bcent,
                        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                        int prob_type, amrex::Box const& vbx)
{
    constexpr amrex::Real pi = 3.1415926535897932;
    amrex::IntVect iv(AMREX_D_DECL(i,j,k));

    if (vbx.contains(iv)) {
        if (flag(i,j,k).isCovered()) {
            phie(i,j,k) = 0.0;
            phib(i,j,k) = 0.0;
        } else {
            amrex::Real x = (i+0.5)*dx[0] - 0.5;
            amrex::Real y = (j+0.5)*dx[1] - 0.5;
            amrex::Real theta = std::atan2(x,y) + 0.5*pi;
            amrex::Real r2 = x*x+y*y;
            phie(i,j,k) = r2*r2 * std::cos(3.*theta);

            x += cent(i,j,k,0)*dx[0];
            y += cent(i,j,k,1)*dx[1];
            theta = std::atan2(x,y) + 0.5*pi;
            r2 = x*x+y*y;
            if (prob_type == 1) {
                rhs(i,j,k) = -7.0 * r2 * std::cos(3.*theta);
            } else {
                rhs(i,j,k) = -(7.0 * r2 - 15. * r2*r2) * std::cos(3.*theta);
            }

            if (flag(i,j,k).isSingleValued()) {
                x = (i+0.5+bcent(i,j,k,0))*dx[0] - 0.5;
                y = (j+0.5+bcent(i,j,k,1))*dx[1] - 0.5;
                theta = std::atan2(x,y) + 0.5*pi;
                r2 = x*x+y*y;
                phib(i,j,k) = r2*r2 * std::cos(3.*theta);
                if (prob_type == 2) {
                    bb(i,j,k) = 1.0 - r2;
                }
            } else {
                phib(i,j,k) = 0.0;
            }
        }
    }

    if (prob_type == 2) {
        if (amrex::surroundingNodes(vbx,0).contains(iv)) {
            amrex::Real x = (i    )*dx[0] - 0.5;
            amrex::Real y = (j+0.5)*dx[1] - 0.5;
            amrex::Real theta = std::atan2(x,y) + 0.5*pi;
            amrex::Real r2 = x*x+y*y;
            bx(i,j,k) = 1.0 - r2;
        }
        if (amrex::surroundingNodes(vbx,1).contains(iv)) {
            amrex::Real x = (i+0.5)*dx[0] - 0.5;
            amrex::Real y = (j    )*dx[1] - 0.5;
            amrex::Real theta = std::atan2(x,y) + 0.5*pi;
            amrex::Real r2 = x*x+y*y;
            by(i,j,k) = 1.0 - r2;
        }
#if (AMREX_SPACEDIM == 3)
        if (amrex::surroundingNodes(vbx,2).contains(iv)) {
            amrex::Real x = (i+0.5)*dx[0] - 0.5;
            amrex::Real y = (j+0.5)*dx[1] - 0.5;
            amrex::Real theta = std::atan2(x,y) + 0.5*pi;
            amrex::Real r2 = x*x+y*y;
            bz(i,j,k) = 1.0 - r2;
        }
#endif
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mytest_set_fluxeb_exact (
    int i, int j, int k,
    amrex::Array4<amrex::Real> const& fluxeb_exact,
    amrex::Array4<amrex::EBCellFlag const> const& flag_arr,
    amrex::EBData const& ebdata,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
    int prob_type, amrex::Box const& vbx)
{
    constexpr amrex::Real pi = 3.1415926535897932;
    amrex::IntVect iv(AMREX_D_DECL(i,j,k));

    if (vbx.contains(iv)) {
        const amrex::EBCellFlag flag = flag_arr(i,j,k);

        if (flag.isCovered()) {
            fluxeb_exact(i,j,k) = 0.0;
        } else if (flag.isSingleValued()) {
            const amrex::Real bcentx = ebdata.get<amrex::EBData_t::bndrycent>(i,j,k,0);
            const amrex::Real bcenty = ebdata.get<amrex::EBData_t::bndrycent>(i,j,k,1);

            // Compute boundary area normal vector
            const amrex::Real apxm = ebdata.get<amrex::EBData_t::apx>(i  ,j,k);
            const amrex::Real apxp = ebdata.get<amrex::EBData_t::apx>(i+1,j,k);
            amrex::Real dapx = (apxm-apxp)/dx[0];
            const amrex::Real apym = ebdata.get<amrex::EBData_t::apy>(i,j  ,k);
            const amrex::Real apyp = ebdata.get<amrex::EBData_t::apy>(i,j+1,k);
            amrex::Real dapy = (apym-apyp)/dx[1];
#if (AMREX_SPACEDIM == 3)
            const amrex::Real apzm = ebdata.get<amrex::EBData_t::apz>(i,j,k  );
            const amrex::Real apzp = ebdata.get<amrex::EBData_t::apz>(i,j,k+1);
            amrex::Real dapz = (apzm-apzp)/dx[2];
#endif
            amrex::Real anorm = std::sqrt(AMREX_D_TERM(dapx*dapx, + dapy*dapy, + dapz*dapz));
            amrex::Real anorminv = amrex::Real(1.0)/anorm;
            amrex::Real nx = dapx * anorminv;
            amrex::Real ny = dapy * anorminv;

            amrex::Real x = (i+0.5+bcentx)*dx[0] - 0.5;
            amrex::Real y = (j+0.5+bcenty)*dx[1] - 0.5;
            amrex::Real theta = std::atan2(x,y) + 0.5*pi;
            amrex::Real r2 = x*x+y*y;

            // Derivation of grad phi Cartesian components:
            // th = atan2(y,x) is the polar angle
            // theta = atan2(x,y) + pi/2
            //       = pi/2 - atan2(y,x) + pi/2
            //       = pi - th
            // phi = r^4 cos(3 theta)
            //     = r^4 cos(3 (pi - th))
            // grad phi = rhat (+4 r^3 cos(3 theta)) + thhat (+3 r^3 sin(3 theta))
            // rhat  = (+x/r) xhat + (+y/r) yhat
            // thhat = (-y/r) xhat + (+x/r) yhat
            // grad_x phi = (+x/r)(+4 r^3 cos(3 theta)) + (-y/r)(+3 r^3 sin(3 theta))
            // grad_y phi = (-y/r)(+4 r^3 cos(3 theta)) + (+x/r)(+3 r^3 sin(3 theta))
            amrex::Real gradphi_r_over_r = 4.0 * r2 * std::cos(3.*theta);
            amrex::Real gradphi_t_over_r = 3.0 * r2 * std::sin(3.*theta);
            amrex::Real gradphi_x = x * gradphi_r_over_r - y * gradphi_t_over_r;
            amrex::Real gradphi_y = y * gradphi_r_over_r + x * gradphi_t_over_r;

            amrex::Real bx = 1.0;
            amrex::Real by = 1.0;
            if (prob_type == 2) {
                bx = 1.0 - r2;
                by = 1.0 - r2;
            }
            // Flux on EB surface = - (area normal).(b tensor).(grad phi)
            amrex::Real n_b_gradphi = nx * bx * gradphi_x + ny * by * gradphi_y;
            fluxeb_exact(i,j,k) = -n_b_gradphi;
        } else {
            fluxeb_exact(i,j,k) = 0.0;
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mytest_set_phi_boundary (int i, int j, int k, amrex::Array4<amrex::Real> const& phi,
                              amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                              amrex::Box const& dbx)
{
    constexpr amrex::Real pi = 3.1415926535897932;

    if (!dbx.contains(amrex::IntVect(AMREX_D_DECL(i,j,k)))) {
        amrex::Real x, y;
        const auto dlo = amrex::lbound(dbx);
        const auto dhi = amrex::ubound(dbx);
        if (i < dlo.x) {
            x = -0.5;
            y = (j+0.5)*dx[1] - 0.5;
        } else if (i > dhi.x) {
            x = 0.5;
            y = (j+0.5)*dx[1] - 0.5;
        } else if (j < dlo.y) {
            x = (i+0.5)*dx[0] - 0.5;
            y = -0.5;

        } else if (j > dhi.y) {
            x = (i+0.5)*dx[0] - 0.5;
            y = 0.5;
        } else {
            x = (i+0.5)*dx[0] - 0.5;
            y = (j+0.5)*dx[1] - 0.5;
        }

        amrex::Real theta = std::atan2(x,y) + 0.5*pi;
        amrex::Real r2 = x*x+y*y;
        phi(i,j,k) = r2*r2 * std::cos(3.*theta);
    }
}

#endif
