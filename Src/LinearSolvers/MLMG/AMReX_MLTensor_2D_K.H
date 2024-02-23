#ifndef AMREX_MLTENSOR_2D_K_H_
#define AMREX_MLTENSOR_2D_K_H_
#include <AMReX_Config.H>

#include <AMReX_MLLinOp_K.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_corners (int icorner, Box const& vbox, // vbox: the valid box
                            Array4<Real> const& vel,
                            Array4<int const> const& mxlo,
                            Array4<int const> const& mylo,
                            Array4<int const> const& mxhi,
                            Array4<int const> const& myhi,
                            Array4<Real const> const& bcvalxlo,
                            Array4<Real const> const& bcvalylo,
                            Array4<Real const> const& bcvalxhi,
                            Array4<Real const> const& bcvalyhi,
                            Array2D<BoundCond,
                                    0,2*AMREX_SPACEDIM,
                                    0,AMREX_SPACEDIM> const& bct,
                            Array2D<Real,
                                    0,2*AMREX_SPACEDIM,
                                    0,AMREX_SPACEDIM> const& bcl,
                            int inhomog, int maxorder,
                            GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                            Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    constexpr int k = 0;
    const auto blen = amrex::length(vbox);
    const auto vlo  = amrex::lbound(vbox);
    const auto vhi  = amrex::ubound(vbox);

    if (icorner == 0) { // xlo & ylo
        int const i = vlo.x-1;
        int const j = vlo.y-1;
        if (mxlo(i,j,k) != BndryData::covered && (dlo.x != vlo.x || dlo.y != vlo.y)) {
            bool x_interior = mylo(i+1,j  ,k) == BndryData::covered; // i+1,j is a valid cell inside domain
            bool x_exterior = mylo(i+1,j  ,k) == BndryData::not_covered; // i+1,j is a ghost cell inside domain
            bool y_interior = mxlo(i  ,j+1,k) == BndryData::covered;
            bool y_exterior = mxlo(i  ,j+1,k) == BndryData::not_covered;
            if ((x_interior && y_interior) || (x_exterior && y_exterior)) {
                for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                }
            } else if (x_interior || dlo.x == vlo.x) {
                for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                }
            } else if (y_interior || dlo.y == vlo.y) {
                for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                }
            }
        }
    } else if (icorner == 1) { // xhi & ylo
        int const i = vhi.x+1;
        int const j = vlo.y-1;
        if (mxhi(i,j,k) != BndryData::covered && (dhi.x != vhi.x || dlo.y != vlo.y)) {
            bool x_interior = mylo(i-1,j  ,k) == BndryData::covered;
            bool x_exterior = mylo(i-1,j  ,k) == BndryData::not_covered;
            bool y_interior = mxhi(i  ,j+1,k) == BndryData::covered;
            bool y_exterior = mxhi(i  ,j+1,k) == BndryData::not_covered;
            if ((x_interior && y_interior) || (x_exterior && y_exterior)) {
                for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                }
            } else if (x_interior || dhi.x == vhi.x) {
                for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                }
            } else if (y_interior || dlo.y == vlo.y) {
                for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                }
            }
        }
    } else if (icorner == 2) { // xlo & yhi
        int const i = vlo.x-1;
        int const j = vhi.y+1;
        if (mxlo(i,j,k) != BndryData::covered && (dlo.x != vlo.x || dhi.y != vhi.y)) {
            bool x_interior = myhi(i+1,j  ,k) == BndryData::covered;
            bool x_exterior = myhi(i+1,j  ,k) == BndryData::not_covered;
            bool y_interior = mxlo(i  ,j-1,k) == BndryData::covered;
            bool y_exterior = mxlo(i  ,j-1,k) == BndryData::not_covered;
            if ((x_interior && y_interior) || (x_exterior && y_exterior)) {
                for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                }
            } else if (x_interior || dlo.x == vlo.x) {
                for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                }
            } else if (y_interior || dhi.y == vhi.y) {
                for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                }
            }
        }
    } else if (icorner == 3) { // xhi & yhi
        int const i = vhi.x+1;
        int const j = vhi.y+1;
        if (mxhi(i,j,k) != BndryData::covered && (dhi.x != vhi.x || dhi.y != vhi.y)) {
            bool x_interior = myhi(i-1,j  ,k) == BndryData::covered;
            bool x_exterior = myhi(i-1,j  ,k) == BndryData::not_covered;
            bool y_interior = mxhi(i  ,j-1,k) == BndryData::covered;
            bool y_exterior = mxhi(i  ,j-1,k) == BndryData::not_covered;
            if ((x_interior && y_interior) || (x_exterior && y_exterior)) {
                for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                }
            } else if (x_interior || dhi.x == vhi.x) {
                for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                }
            } else if (y_interior || dhi.y == vhi.y) {
                for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms_fx (Box const& box, Array4<Real> const& fx,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& etax,
                              Array4<Real const> const& kapx,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = Real(2./3.);

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real dudy = mltensor_dy_on_xface(i,j,k,0,vel,dyi);
            Real dvdy = mltensor_dy_on_xface(i,j,k,1,vel,dyi);
            Real divu = dvdy;
            Real xif = kapx(i,j,0);
            Real mun = Real(0.75)*(etax(i,j,0,0)-xif);  // restore the original eta
            Real mut =             etax(i,j,0,1);
            fx(i,j,0,0) = -mun*(-twoThirds*divu) - xif*divu;
            fx(i,j,0,1) = -mut*dudy;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms_fy (Box const& box, Array4<Real> const& fy,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& etay,
                              Array4<Real const> const& kapy,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = Real(2./3.);

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real dudx = mltensor_dx_on_yface(i,j,k,0,vel,dxi);
            Real dvdx = mltensor_dx_on_yface(i,j,k,1,vel,dxi);
            Real divu = dudx;
            Real xif = kapy(i,j,0);
            Real mun = Real(0.75)*(etay(i,j,0,1)-xif);  // restore the original eta
            Real mut =             etay(i,j,0,0);
            fy(i,j,0,0) = -mut*dvdx;
            fy(i,j,0,1) = -mun*(-twoThirds*divu) - xif*divu;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms_fx (Box const& box, Array4<Real> const& fx,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& etax,
                              Array4<Real const> const& kapx,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Array4<Real const> const& bvxlo,
                              Array4<Real const> const& bvxhi,
                              Array2D<BoundCond,
                                      0,2*AMREX_SPACEDIM,
                                      0,AMREX_SPACEDIM> const& bct,
                              Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = Real(2./3.);

    // Three BC types: reflect odd, neumann, and dirichlet

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            Real dudy = mltensor_dy_on_xface(i,j,k,0,vel,dyi,bvxlo,bvxhi,bct,dlo,dhi);
            Real dvdy = mltensor_dy_on_xface(i,j,k,1,vel,dyi,bvxlo,bvxhi,bct,dlo,dhi);
            Real divu = dvdy;
            Real xif = kapx(i,j,0);
            Real mun = Real(0.75)*(etax(i,j,0,0)-xif);  // restore the original eta
            Real mut =             etax(i,j,0,1);
            fx(i,j,0,0) = -mun*(-twoThirds*divu) - xif*divu;
            fx(i,j,0,1) = -mut*dudy;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms_fy (Box const& box, Array4<Real> const& fy,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& etay,
                              Array4<Real const> const& kapy,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Array4<Real const> const& bvylo,
                              Array4<Real const> const& bvyhi,
                              Array2D<BoundCond,
                                      0,2*AMREX_SPACEDIM,
                                      0,AMREX_SPACEDIM> const& bct,
                              Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = Real(2./3.);

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            Real dudx = mltensor_dx_on_yface(i,j,k,0,vel,dxi,bvylo,bvyhi,bct,dlo,dhi);
            Real dvdx = mltensor_dx_on_yface(i,j,k,1,vel,dxi,bvylo,bvyhi,bct,dlo,dhi);
            Real divu = dudx;
            Real xif = kapy(i,j,0);
            Real mun = Real(0.75)*(etay(i,j,0,1)-xif);  // restore the original eta
            Real mut =             etay(i,j,0,0);
            fy(i,j,0,0) = -mut*dvdx;
            fy(i,j,0,1) = -mun*(-twoThirds*divu) - xif*divu;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms (Box const& box, Array4<Real> const& Ax,
                           Array4<Real const> const& fx,
                           Array4<Real const> const& fy,
                           GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                           Real bscalar) noexcept
{
    const Real dxi = bscalar * dxinv[0];
    const Real dyi = bscalar * dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Ax(i,j,0,0) += dxi*(fx(i+1,j  ,0,0) - fx(i,j,0,0))
                +          dyi*(fy(i  ,j+1,0,0) - fy(i,j,0,0));
            Ax(i,j,0,1) += dxi*(fx(i+1,j  ,0,1) - fx(i,j,0,1))
                +          dyi*(fy(i  ,j+1,0,1) - fy(i,j,0,1));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms_os (Box const& box, Array4<Real> const& Ax,
                              Array4<Real const> const& fx,
                              Array4<Real const> const& fy,
                              Array4<int const> const& osm,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Real bscalar) noexcept
{
    const Real dxi = bscalar * dxinv[0];
    const Real dyi = bscalar * dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            if (osm(i,j,0) == 0) {
                Ax(i,j,0,0) = 0.0;
                Ax(i,j,0,1) = 0.0;
            } else {
                Ax(i,j,0,0) += dxi*(fx(i+1,j  ,0,0) - fx(i,j,0,0))
                    +          dyi*(fy(i  ,j+1,0,0) - fy(i,j,0,0));
                Ax(i,j,0,1) += dxi*(fx(i+1,j  ,0,1) - fx(i,j,0,1))
                    +          dyi*(fy(i  ,j+1,0,1) - fy(i,j,0,1));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_vel_grads_fx (Box const& box, Array4<Real> const& fx,
                              Array4<Real const> const& vel,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real dudx = (vel(i,j,0,0) - vel(i-1,j,0,0))*dxi;
            Real dvdx = (vel(i,j,0,1) - vel(i-1,j,0,1))*dxi;
            Real dudy = mltensor_dy_on_xface(i,j,k,0,vel,dyi);
            Real dvdy = mltensor_dy_on_xface(i,j,k,1,vel,dyi);
            fx(i,j,0,0) = dudx;
            fx(i,j,0,1) = dvdx;
            fx(i,j,0,2) = dudy;
            fx(i,j,0,3) = dvdy;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_vel_grads_fy (Box const& box, Array4<Real> const& fy,
                              Array4<Real const> const& vel,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real dudx = mltensor_dx_on_yface(i,j,k,0,vel,dxi);
            Real dvdx = mltensor_dx_on_yface(i,j,k,1,vel,dxi);
            Real dudy = (vel(i,j,0,0) - vel(i,j-1,0,0))*dyi;
            Real dvdy = (vel(i,j,0,1) - vel(i,j-1,0,1))*dyi;
            fy(i,j,0,0) = dudx;
            fy(i,j,0,1) = dvdx;
            fy(i,j,0,2) = dudy;
            fy(i,j,0,3) = dvdy;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_vel_grads_fx (Box const& box, Array4<Real> const& fx,
                            Array4<Real const> const& vel,
                            GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                            Array4<Real const> const& bvxlo,
                            Array4<Real const> const& bvxhi,
                            Array2D<BoundCond,
                                    0,2*AMREX_SPACEDIM,
                                    0,AMREX_SPACEDIM> const& bct,
                            Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            Real dudx = (vel(i,j,0,0) - vel(i-1,j,0,0))*dxi;
            Real dvdx = (vel(i,j,0,1) - vel(i-1,j,0,1))*dxi;
            Real dudy = mltensor_dy_on_xface(i,j,k,0,vel,dyi,bvxlo,bvxhi,bct,dlo,dhi);
            Real dvdy = mltensor_dy_on_xface(i,j,k,1,vel,dyi,bvxlo,bvxhi,bct,dlo,dhi);
            fx(i,j,0,0) = dudx;
            fx(i,j,0,1) = dvdx;
            fx(i,j,0,2) = dudy;
            fx(i,j,0,3) = dvdy;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_vel_grads_fy (Box const& box, Array4<Real> const& fy,
                            Array4<Real const> const& vel,
                            GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                            Array4<Real const> const& bvylo,
                            Array4<Real const> const& bvyhi,
                            Array2D<BoundCond,
                                    0,2*AMREX_SPACEDIM,
                                    0,AMREX_SPACEDIM> const& bct,
                            Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            Real dudx = mltensor_dx_on_yface(i,j,k,0,vel,dxi,bvylo,bvyhi,bct,dlo,dhi);
            Real dvdx = mltensor_dx_on_yface(i,j,k,1,vel,dxi,bvylo,bvyhi,bct,dlo,dhi);
            Real dudy = (vel(i,j,0,0) - vel(i,j-1,0,0))*dyi;
            Real dvdy = (vel(i,j,0,1) - vel(i,j-1,0,1))*dyi;
            fy(i,j,0,0) = dudx;
            fy(i,j,0,1) = dvdx;
            fy(i,j,0,2) = dudy;
            fy(i,j,0,3) = dvdy;
        }
    }
}

}

#endif
