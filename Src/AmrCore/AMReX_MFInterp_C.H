#ifndef AMREX_MF_INTERP_C_H_
#define AMREX_MF_INTERP_C_H_
#include <AMReX_Config.H>

#include <AMReX_Array4.H>
#include <AMReX_BCRec.H>

namespace amrex {

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mf_compute_slopes_x (int i, int j, int k, Array4<Real const> const& u, int nu,
                          Box const& domain, BCRec const& bc)
{
    Real dc = Real(0.5) * (u(i+1,j,k,nu) - u(i-1,j,k,nu));
    if (i == domain.smallEnd(0) && (bc.lo(0) == BCType::ext_dir ||
                                    bc.lo(0) == BCType::hoextrap)) {
        if (i+2 < u.end.x) {
            dc = -Real(16./15.)*u(i-1,j,k,nu) + Real(0.5)*u(i,j,k,nu)
                + Real(2./3.)*u(i+1,j,k,nu) - Real(0.1)*u(i+2,j,k,nu);
        } else {
            dc = Real(0.25)*(u(i+1,j,k,nu)+Real(5.)*u(i,j,k,nu)-Real(6.)*u(i-1,j,k,nu));
        }
    }
    if (i == domain.bigEnd(0) && (bc.hi(0) == BCType::ext_dir ||
                                  bc.hi(0) == BCType::hoextrap)) {
        if (i-2 >= u.begin.x) {
            dc = Real(16./15.)*u(i+1,j,k,nu) - Real(0.5)*u(i,j,k,nu)
                - Real(2./3.)*u(i-1,j,k,nu) + Real(0.1)*u(i-2,j,k,nu);
        } else {
            dc = -Real(0.25)*(u(i-1,j,k,nu)+Real(5.)*u(i,j,k,nu)-Real(6.)*u(i+1,j,k,nu));
        }
    }
    return dc;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mf_compute_slopes_y (int i, int j, int k, Array4<Real const> const& u, int nu,
                          Box const& domain, BCRec const& bc)
{
    Real dc = Real(0.5) * (u(i,j+1,k,nu) - u(i,j-1,k,nu));
    if (j == domain.smallEnd(1) && (bc.lo(1) == BCType::ext_dir ||
                                    bc.lo(1) == BCType::hoextrap)) {
        if (j+2 < u.end.y) {
            dc = -Real(16./15.)*u(i,j-1,k,nu) + Real(0.5)*u(i,j,k,nu)
                + Real(2./3.)*u(i,j+1,k,nu) - Real(0.1)*u(i,j+2,k,nu);
        } else {
            dc = Real(0.25)*(u(i,j+1,k,nu)+Real(5.)*u(i,j,k,nu)-Real(6.)*u(i,j-1,k,nu));
        }
    }
    if (j == domain.bigEnd(1) && (bc.hi(1) == BCType::ext_dir ||
                                  bc.hi(1) == BCType::hoextrap)) {
        if (j-2 >= u.begin.y) {
            dc = Real(16./15.)*u(i,j+1,k,nu) - Real(0.5)*u(i,j,k,nu)
                - Real(2./3.)*u(i,j-1,k,nu) + Real(0.1)*u(i,j-2,k,nu);
        } else {
            dc = -Real(0.25)*(u(i,j-1,k,nu)+Real(5.)*u(i,j,k,nu)-Real(6.)*u(i,j+1,k,nu));
        }
    }
    return dc;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mf_compute_slopes_z (int i, int j, int k, Array4<Real const> const& u, int nu,
                          Box const& domain, BCRec const& bc)
{
    Real dc = Real(0.5) * (u(i,j,k+1,nu) - u(i,j,k-1,nu));
    if (k == domain.smallEnd(2) && (bc.lo(2) == BCType::ext_dir ||
                                    bc.lo(2) == BCType::hoextrap)) {
        if (k+2 < u.end.z) {
            dc = -Real(16./15.)*u(i,j,k-1,nu) + Real(0.5)*u(i,j,k,nu)
                + Real(2./3.)*u(i,j,k+1,nu) - Real(0.1)*u(i,j,k+2,nu);
        } else {
            dc = Real(0.25)*(u(i,j,k+1,nu)+Real(5.)*u(i,j,k,nu)-Real(6.)*u(i,j,k-1,nu));
        }
    }
    if (k == domain.bigEnd(2) && (bc.hi(2) == BCType::ext_dir ||
                                  bc.hi(2) == BCType::hoextrap)) {
        if (k-2 >= u.begin.z) {
            dc = Real(16./15.)*u(i,j,k+1,nu) - Real(0.5)*u(i,j,k,nu)
                - Real(2./3.)*u(i,j,k-1,nu) + Real(0.1)*u(i,j,k-2,nu);
        } else {
            dc = -Real(0.25)*(u(i,j,k-1,nu)+Real(5.)*u(i,j,k,nu)-Real(6.)*u(i,j,k+1,nu));
        }
    }
    return dc;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mf_cell_quadratic_compute_slopes_xx (int i, int j, int k,
                                          Array4<Real const> const& u, int nu,
                                          Box const& domain, BCRec const& bc)
{
    Real xx = u(i-1,j,k,nu) - 2.0_rt * u(i,j,k,nu) + u(i+1,j,k,nu);
    if (i == domain.smallEnd(0) && (bc.lo(0) == BCType::ext_dir ||
                                    bc.lo(0) == BCType::hoextrap)) {
        if (i+2 < u.end.x) {
            xx = 0._rt;
        }
    }
    if (i == domain.bigEnd(0) && (bc.hi(0) == BCType::ext_dir ||
                                  bc.hi(0) == BCType::hoextrap)) {
        if (i-2 >= u.begin.x) {
            xx = 0._rt;
        }
    }
    return xx;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mf_cell_quadratic_compute_slopes_yy (int i, int j, int k,
                                          Array4<Real const> const& u, int nu,
                                          Box const& domain, BCRec const& bc)
{
    Real yy = u(i,j-1,k,nu) - 2.0_rt * u(i,j,k,nu) + u(i,j+1,k,nu);
    if (j == domain.smallEnd(1) && (bc.lo(1) == BCType::ext_dir ||
                                    bc.lo(1) == BCType::hoextrap)) {
        if (j+2 < u.end.y) {
            yy = 0._rt;
        }
    }
    if (j == domain.bigEnd(1) && (bc.hi(1) == BCType::ext_dir ||
                                  bc.hi(1) == BCType::hoextrap)) {
        if (j-2 >= u.begin.y) {
            yy = 0._rt;
        }
    }
    return yy;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mf_cell_quadratic_compute_slopes_zz (int i, int j, int k,
                                          Array4<Real const> const& u, int nu,
                                          Box const& domain, BCRec const& bc)
{
    Real zz = u(i,j,k-1,nu) - 2.0_rt * u(i,j,k,nu) + u(i,j,k+1,nu);
    if (k == domain.smallEnd(2) && (bc.lo(2) == BCType::ext_dir ||
                                    bc.lo(2) == BCType::hoextrap)) {
        if (k+2 < u.end.z) {
            zz = 0._rt;
        }
    }
    if (k == domain.bigEnd(1) && (bc.hi(2) == BCType::ext_dir ||
                                  bc.hi(2) == BCType::hoextrap)) {
        if (k-2 >= u.begin.z) {
            zz = 0._rt;
        }
    }
    return zz;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mf_cell_quadratic_compute_slopes_xy (int i, int j, int k,
                                          Array4<Real const> const& u, int nu,
                                          Box const& domain, BCRec const& bc)
{
    Real xy = (1._rt/4._rt) * ( u(i-1,j-1,k,nu) - u(i+1,j-1,k,nu)
                              - u(i-1,j+1,k,nu) + u(i+1,j+1,k,nu) );
    if (i == domain.smallEnd(0) && (bc.lo(0) == BCType::ext_dir ||
                                    bc.lo(0) == BCType::hoextrap)) {
        if (i+2 < u.end.x) {
            xy = 0._rt;
        }
    }
    if (i == domain.bigEnd(0) && (bc.hi(0) == BCType::ext_dir ||
                                  bc.hi(0) == BCType::hoextrap)) {
        if (i-2 >= u.begin.x) {
            xy = 0._rt;
        }
    }
    if (j == domain.smallEnd(1) && (bc.lo(1) == BCType::ext_dir ||
                                    bc.lo(1) == BCType::hoextrap)) {
        if (j+2 < u.end.y) {
            xy = 0._rt;
        }
    }
    if (j == domain.bigEnd(1) && (bc.hi(1) == BCType::ext_dir ||
                                  bc.hi(1) == BCType::hoextrap)) {
        if (j-2 >= u.begin.y) {
            xy = 0._rt;
        }
    }
    return xy;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mf_cell_quadratic_compute_slopes_xz (int i, int j, int k,
                                          Array4<Real const> const& u, int nu,
                                          Box const& domain, BCRec const& bc)
{
    Real xz = (1._rt/4._rt) * ( u(i-1,j,k-1,nu) - u(i+1,j,k-1,nu)
                              - u(i-1,j,k+1,nu) + u(i+1,j,k+1,nu) );
    if (i == domain.smallEnd(0) && (bc.lo(0) == BCType::ext_dir ||
                                    bc.lo(0) == BCType::hoextrap)) {
        if (i+2 < u.end.x) {
            xz = 0._rt;
        }
    }
    if (i == domain.bigEnd(0) && (bc.hi(0) == BCType::ext_dir ||
                                  bc.hi(0) == BCType::hoextrap)) {
        if (i-2 >= u.begin.x) {
            xz = 0._rt;
        }
    }
    if (k == domain.smallEnd(2) && (bc.lo(2) == BCType::ext_dir ||
                                    bc.lo(2) == BCType::hoextrap)) {
        if (k+2 < u.end.z) {
            xz = 0._rt;
        }
    }
    if (k == domain.bigEnd(1) && (bc.hi(2) == BCType::ext_dir ||
                                  bc.hi(2) == BCType::hoextrap)) {
        if (k-2 >= u.begin.z) {
            xz = 0._rt;
        }
    }
    return xz;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mf_cell_quadratic_compute_slopes_yz (int i, int j, int k,
                                          Array4<Real const> const& u, int nu,
                                          Box const& domain, BCRec const& bc)
{
    Real yz = (1._rt/4._rt) * ( u(i,j-1,k-1,nu) - u(i,j-1,k+1,nu)
                              - u(i,j+1,k-1,nu) + u(i,j+1,k+1,nu) );
    if (j == domain.smallEnd(1) && (bc.lo(1) == BCType::ext_dir ||
                                    bc.lo(1) == BCType::hoextrap)) {
        if (j+2 < u.end.y) {
            yz = 0._rt;
        }
    }
    if (j == domain.bigEnd(1) && (bc.hi(1) == BCType::ext_dir ||
                                  bc.hi(1) == BCType::hoextrap)) {
        if (j-2 >= u.begin.y) {
            yz = 0._rt;
        }
    }
    if (k == domain.smallEnd(2) && (bc.lo(2) == BCType::ext_dir ||
                                    bc.lo(2) == BCType::hoextrap)) {
        if (k+2 < u.end.z) {
            yz = 0._rt;
        }
    }
    if (k == domain.bigEnd(1) && (bc.hi(2) == BCType::ext_dir ||
                                  bc.hi(2) == BCType::hoextrap)) {
        if (k-2 >= u.begin.z) {
            yz = 0._rt;
        }
    }
    return yz;
}

} // namespace

} // namespace amrex

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MFInterp_1D_C.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MFInterp_2D_C.H>
#else
#include <AMReX_MFInterp_3D_C.H>
#endif

#endif
