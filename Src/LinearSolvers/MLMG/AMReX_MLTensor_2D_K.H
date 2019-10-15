#ifndef AMREX_MLTENSOR_2D_K_H_
#define AMREX_MLTENSOR_2D_K_H_

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
                            GpuArray<BoundCond,2*AMREX_SPACEDIM*AMREX_SPACEDIM> const& bct,
                            GpuArray<Real,2*AMREX_SPACEDIM*AMREX_SPACEDIM> const& bcl,
                            int inhomog, int maxorder,
                            GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                            Box const& domain) noexcept
{
    constexpr int oxlo = 0;
    constexpr int oylo = 1;
    constexpr int oxhi = 2;
    constexpr int oyhi = 3;
    constexpr int xdir = 0;
    constexpr int ydir = 1;
    const auto blen = amrex::length(vbox);
    const auto vlo  = amrex::lbound(vbox);
    const auto vhi  = amrex::ubound(vbox);
    const auto dlo  = amrex::lbound(domain);
    const auto dhi  = amrex::ubound(domain);
    for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
        switch (icorner) {
        case 0: {
            // xlo & ylo
            Box bx = amrex::adjCellLo(amrex::adjCellLo(vbox,xdir,1),ydir,1);
            if (vlo.x == dlo.x and vlo.y == dlo.y) {
                vel(vlo.x-1,vlo.y-1,0,icomp) = vel(vlo.x-1,vlo.y,0,icomp)
                    + vel(vlo.x,vlo.y-1,0,icomp) - vel(vlo.x,vlo.y,0,icomp);
            } else if (vlo.x == dlo.x) {
                int offset = AMREX_SPACEDIM * oxlo;
                mllinop_apply_bc_x(Orientation::low, bx, blen.x,
                                   vel, mxlo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxlo, maxorder, dxinv[xdir], inhomog, icomp);
            } else if (vlo.y == dlo.y) {
                int offset = AMREX_SPACEDIM * oylo;
                mllinop_apply_bc_y(Orientation::low, bx, blen.y,
                                   vel, mylo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalylo, maxorder, dxinv[ydir], inhomog, icomp);
            }
            break;
        }
        case 1: {
            // xhi & ylo
            Box bx = amrex::adjCellLo(amrex::adjCellHi(vbox,xdir,1),ydir,1);
            if (vhi.x == dhi.x and vlo.y == dlo.y) {
                vel(vhi.x+1,vlo.y-1,0,icomp) = vel(vhi.x+1,vlo.y,0,icomp)
                    + vel(vhi.x,vlo.y-1,0,icomp) - vel(vhi.x,vlo.y,0,icomp);
            } else if (vhi.x == dhi.x) {
                int offset = AMREX_SPACEDIM * oxhi;
                mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                   vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
            } else if (vlo.y == dlo.y) {
                int offset = AMREX_SPACEDIM * oylo;
                mllinop_apply_bc_y(Orientation::low, bx, blen.y,
                                   vel, mylo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalylo, maxorder, dxinv[ydir], inhomog, icomp);
            }
            break;
        }
        case 2: {
            // xlo & yhi
            Box bx = amrex::adjCellHi(amrex::adjCellLo(vbox,xdir,1),ydir,1);
            if (vlo.x == dlo.x and vhi.y == dhi.y) {
                vel(vlo.x-1,vhi.y+1,0,icomp) = vel(vlo.x-1,vhi.y,0,icomp)
                    + vel(vlo.x,vhi.y+1,0,icomp) - vel(vlo.x,vhi.y,0,icomp);
            } else if (vlo.x == dlo.x) {
                int offset = AMREX_SPACEDIM * oxlo;
                mllinop_apply_bc_x(Orientation::low, bx, blen.x,
                                   vel, mxlo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxlo, maxorder, dxinv[xdir], inhomog, icomp);
            } else if (vhi.y == dhi.y) {
                int offset = AMREX_SPACEDIM * oyhi;
                mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                   vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
            }
            break;
        }
        case 3: {
            // xhi & yhi
            Box bx = amrex::adjCellHi(amrex::adjCellHi(vbox,xdir,1),ydir,1);
            if (vhi.x == dhi.x and vhi.y == dhi.y) {
                vel(vhi.x+1,vhi.y+1,0,icomp) = vel(vhi.x+1,vhi.y,0,icomp)
                    + vel(vhi.x,vhi.y+1,0,icomp) - vel(vhi.x,vhi.y,0,icomp);
            } else if (vhi.x == dhi.x) {
                int offset = AMREX_SPACEDIM * oxhi;
                mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                   vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
            } else  if (vhi.y == dhi.y) {
                int offset = AMREX_SPACEDIM * oyhi;
                mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                   vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
            }
            break;
        }
        default: {}
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
    constexpr Real twoThirds = 2./3.;

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real dudy = (vel(i,j+1,0,0)+vel(i-1,j+1,0,0)-vel(i,j-1,0,0)-vel(i-1,j-1,0,0))*(0.25*dyi);
            Real dvdy = (vel(i,j+1,0,1)+vel(i-1,j+1,0,1)-vel(i,j-1,0,1)-vel(i-1,j-1,0,1))*(0.25*dyi);
            Real divu = dvdy;
            Real xif = kapx(i,j,0);
            Real mun = 0.75*(etax(i,j,0,0)-xif);  // restore the original eta
            Real mut =       etax(i,j,0,1);
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
    constexpr Real twoThirds = 2./3.;

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real dudx = (vel(i+1,j,0,0)+vel(i+1,j-1,0,0)-vel(i-1,j,0,0)-vel(i-1,j-1,0,0))*(0.25*dxi);
            Real dvdx = (vel(i+1,j,0,1)+vel(i+1,j-1,0,1)-vel(i-1,j,0,1)-vel(i-1,j-1,0,1))*(0.25*dxi);
            Real divu = dudx;
            Real xif = kapy(i,j,0);
            Real mun = 0.75*(etay(i,j,0,1)-xif);  // restore the original eta
            Real mut =       etay(i,j,0,0);
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

}

#endif
