#ifndef AMREX_MLTENSOR_3D_K_H_
#define AMREX_MLTENSOR_3D_K_H_

#include <AMReX_MLLinOp_K.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_corners (int icorner, Box const& vbox, // vbox: the valid box
                            Array4<Real> const& vel,
                            Array4<int const> const& mxlo,
                            Array4<int const> const& mylo,
                            Array4<int const> const& mzlo,
                            Array4<int const> const& mxhi,
                            Array4<int const> const& myhi,
                            Array4<int const> const& mzhi,
                            Array4<Real const> const& bcvalxlo,
                            Array4<Real const> const& bcvalylo,
                            Array4<Real const> const& bcvalzlo,
                            Array4<Real const> const& bcvalxhi,
                            Array4<Real const> const& bcvalyhi,
                            Array4<Real const> const& bcvalzhi,
                            GpuArray<BoundCond,2*AMREX_SPACEDIM*AMREX_SPACEDIM> const& bct,
                            GpuArray<Real,2*AMREX_SPACEDIM*AMREX_SPACEDIM> const& bcl,
                            int inhomog, int maxorder,
                            GpuArray<Real,AMREX_SPACEDIM> const& dxinv, Box const& domain) noexcept
{
    constexpr int oxlo = 0;
    constexpr int oylo = 1;
    constexpr int ozlo = 2;
    constexpr int oxhi = 3;
    constexpr int oyhi = 4;
    constexpr int ozhi = 5;
    constexpr int xdir = 0;
    constexpr int ydir = 1;
    constexpr int zdir = 2;
    const auto blen = amrex::length(vbox);
    const auto vlo  = amrex::lbound(vbox);
    const auto vhi  = amrex::ubound(vbox);
    const auto dlo  = amrex::lbound(domain);
    const auto dhi  = amrex::ubound(domain);
    for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
        switch (icorner) {
        case 0: {
            // xlo & ylo & zlo
            Box bx = amrex::adjCellLo(amrex::adjCellLo(amrex::adjCellLo(vbox,xdir,1),ydir,1),zdir,1);
            if (vlo.x == dlo.x and vlo.y == dlo.y and vlo.z == dlo.z) {
                vel      (vlo.x-1,vlo.y-1,vlo.z-1,icomp)
                    = vel(vlo.x-1,vlo.y  ,vlo.z  ,icomp)
                    + vel(vlo.x  ,vlo.y-1,vlo.z  ,icomp)
                    + vel(vlo.x  ,vlo.y  ,vlo.z-1,icomp)
                    - vel(vlo.x  ,vlo.y  ,vlo.z  ,icomp) * 2.0;
            } else if (vlo.x == dlo.x and vlo.y == dlo.y) {
                vel      (vlo.x-1,vlo.y-1,vlo.z-1,icomp)
                    = vel(vlo.x-1,vlo.y  ,vlo.z-1,icomp)
                    + vel(vlo.x  ,vlo.y-1,vlo.z-1,icomp)
                    - vel(vlo.x  ,vlo.y  ,vlo.z-1,icomp);
            } else if (vlo.x == dlo.x and vlo.z == dlo.z) {
                vel      (vlo.x-1,vlo.y-1,vlo.z-1,icomp)
                    = vel(vlo.x-1,vlo.y-1,vlo.z  ,icomp)
                    + vel(vlo.x  ,vlo.y-1,vlo.z-1,icomp)
                    - vel(vlo.x  ,vlo.y-1,vlo.z  ,icomp);
            } else if (vlo.y == dlo.y and vlo.z == dlo.z) {
                vel      (vlo.x-1,vlo.y-1,vlo.z-1,icomp)
                    = vel(vlo.x-1,vlo.y-1,vlo.z  ,icomp)
                    + vel(vlo.x-1,vlo.y  ,vlo.z-1,icomp)
                    - vel(vlo.x-1,vlo.y  ,vlo.z  ,icomp);
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
            } else if (vlo.z == dlo.z) {
                int offset = AMREX_SPACEDIM * ozlo;
                mllinop_apply_bc_z(Orientation::low, bx, blen.z,
                                   vel, mzlo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalzlo, maxorder, dxinv[zdir], inhomog, icomp);
            } else if (mxlo(vlo.x-1,vlo.y-1,vlo.z-1) != BndryData::covered) {
                if (mylo(vlo.x,vlo.y-1,vlo.z-1) == BndryData::covered) {
                    int offset = AMREX_SPACEDIM * oxlo;
                    mllinop_apply_bc_x(Orientation::low, bx, blen.x,
                                       vel, mxlo, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalxlo, maxorder, dxinv[xdir], inhomog, icomp);
                } else if (mxlo(vlo.x-1,vlo.y,vlo.z-1) == BndryData::covered) {
                    int offset = AMREX_SPACEDIM * oylo;
                    mllinop_apply_bc_y(Orientation::low, bx, blen.y,
                                       vel, mylo, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalylo, maxorder, dxinv[ydir], inhomog, icomp);
                } else {
                    int offset = AMREX_SPACEDIM * ozlo;
                    mllinop_apply_bc_z(Orientation::low, bx, blen.z,
                                       vel, mzlo, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalzlo, maxorder, dxinv[zdir], inhomog, icomp);
                }
            }
            break;
        }
        case 1: {
            // xhi & ylo & zlo
            Box bx = amrex::adjCellLo(amrex::adjCellLo(amrex::adjCellHi(vbox,xdir,1),ydir,1),zdir,1);
            if (vhi.x == dhi.x and vlo.y == dlo.y and vlo.z == dlo.z) {
                vel      (vhi.x+1,vlo.y-1,vlo.z-1,icomp)
                    = vel(vhi.x+1,vlo.y  ,vlo.z  ,icomp)
                    + vel(vhi.x  ,vlo.y-1,vlo.z  ,icomp)
                    + vel(vhi.x  ,vlo.y  ,vlo.z-1,icomp)
                    - vel(vhi.x  ,vlo.y  ,vlo.z  ,icomp) * 2.0;
            } else if (vhi.x == dhi.x and vlo.y == dlo.y) {
                vel      (vhi.x+1,vlo.y-1,vlo.z-1,icomp)
                    = vel(vhi.x+1,vlo.y  ,vlo.z-1,icomp)
                    + vel(vhi.x  ,vlo.y-1,vlo.z-1,icomp)
                    - vel(vhi.x  ,vlo.y  ,vlo.z-1,icomp);
            } else if (vhi.x == dhi.x and vlo.z == dlo.z) {
                vel      (vhi.x+1,vlo.y-1,vlo.z-1,icomp)
                    = vel(vhi.x+1,vlo.y-1,vlo.z  ,icomp)
                    + vel(vhi.x  ,vlo.y-1,vlo.z-1,icomp)
                    - vel(vhi.x  ,vlo.y-1,vlo.z  ,icomp);
            } else if (vlo.y == dlo.y and vlo.z == dlo.z) {
                vel      (vhi.x+1,vlo.y-1,vlo.z-1,icomp)
                    = vel(vhi.x+1,vlo.y-1,vlo.z  ,icomp)
                    + vel(vhi.x+1,vlo.y  ,vlo.z-1,icomp)
                    - vel(vhi.x+1,vlo.y  ,vlo.z  ,icomp);
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
            } else if (vlo.z == dlo.z) {
                int offset = AMREX_SPACEDIM * ozlo;
                mllinop_apply_bc_z(Orientation::low, bx, blen.z,
                                   vel, mzlo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalzlo, maxorder, dxinv[zdir], inhomog, icomp);
            } else if (mxhi(vhi.x+1,vlo.y-1,vlo.z-1) != BndryData::covered) {
                if (mylo(vhi.x,vlo.y-1,vlo.z-1) == BndryData::covered) {
                    int offset = AMREX_SPACEDIM * oxhi;
                    mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                       vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
                } else if (mxhi(vhi.x+1,vlo.y,vlo.z-1) == BndryData::covered) {
                    int offset = AMREX_SPACEDIM * oylo;
                    mllinop_apply_bc_y(Orientation::low, bx, blen.y,
                                       vel, mylo, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalylo, maxorder, dxinv[ydir], inhomog, icomp);
                } else {
                    int offset = AMREX_SPACEDIM * ozlo;
                    mllinop_apply_bc_z(Orientation::low, bx, blen.z,
                                       vel, mzlo, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalzlo, maxorder, dxinv[zdir], inhomog, icomp);
                }
            }
            break;
        }
        case 2: {
            // xlo & yhi & zlo
            Box bx = amrex::adjCellLo(amrex::adjCellHi(amrex::adjCellLo(vbox,xdir,1),ydir,1),zdir,1);
            if (vlo.x == dlo.x and vhi.y == dhi.y and vlo.z == dlo.z) {
                vel      (vlo.x-1,vhi.y+1,vlo.z-1,icomp)
                    = vel(vlo.x-1,vhi.y  ,vlo.z  ,icomp)
                    + vel(vlo.x  ,vhi.y+1,vlo.z  ,icomp)
                    + vel(vlo.x  ,vhi.y  ,vlo.z-1,icomp)
                    - vel(vlo.x  ,vhi.y  ,vlo.z  ,icomp) * 2.0;
            } else if (vlo.x == dlo.x and vhi.y == dhi.y) {
                vel      (vlo.x-1,vhi.y+1,vlo.z-1,icomp)
                    = vel(vlo.x-1,vhi.y  ,vlo.z-1,icomp)
                    + vel(vlo.x  ,vhi.y+1,vlo.z-1,icomp)
                    - vel(vlo.x  ,vhi.y  ,vlo.z-1,icomp);
            } else if (vlo.x == dlo.x and vlo.z == dlo.z) {
                vel      (vlo.x-1,vhi.y+1,vlo.z-1,icomp)
                    = vel(vlo.x-1,vhi.y+1,vlo.z  ,icomp)
                    + vel(vlo.x  ,vhi.y+1,vlo.z-1,icomp)
                    - vel(vlo.x  ,vhi.y+1,vlo.z  ,icomp);
            } else if (vhi.y == dhi.y and vlo.z == dlo.z) {
                vel      (vlo.x-1,vhi.y+1,vlo.z-1,icomp)
                    = vel(vlo.x-1,vhi.y+1,vlo.z  ,icomp)
                    + vel(vlo.x-1,vhi.y  ,vlo.z-1,icomp)
                    - vel(vlo.x-1,vhi.y  ,vlo.z  ,icomp);
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
            } else if (vlo.z == dlo.z) {
                int offset = AMREX_SPACEDIM * ozlo;
                mllinop_apply_bc_z(Orientation::low, bx, blen.z,
                                   vel, mzlo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalzlo, maxorder, dxinv[zdir], inhomog, icomp);
            } else if (mxlo(vlo.x-1,vhi.y+1,vlo.z-1) != BndryData::covered) {
                if (myhi(vlo.x,vhi.y+1,vlo.z-1) == BndryData::covered) {
                    int offset = AMREX_SPACEDIM * oxlo;
                    mllinop_apply_bc_x(Orientation::low, bx, blen.x,
                                       vel, mxlo, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalxlo, maxorder, dxinv[xdir], inhomog, icomp);
                } else if (mxlo(vlo.x-1,vhi.y,vlo.z-1) == BndryData::covered) {
                    int offset = AMREX_SPACEDIM * oyhi;
                    mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                       vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
                } else {
                    int offset = AMREX_SPACEDIM * ozlo;
                    mllinop_apply_bc_z(Orientation::low, bx, blen.z,
                                       vel, mzlo, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalzlo, maxorder, dxinv[zdir], inhomog, icomp);
                }
            }
            break;
        }
        case 3: {
            // xhi & yhi & zlo
            Box bx = amrex::adjCellLo(amrex::adjCellHi(amrex::adjCellHi(vbox,xdir,1),ydir,1),zdir,1);
            if (vhi.x == dhi.x and vhi.y == dhi.y and vlo.z == dlo.z) {
                vel      (vhi.x+1,vhi.y+1,vlo.z-1,icomp)
                    = vel(vhi.x+1,vhi.y  ,vlo.z  ,icomp)
                    + vel(vhi.x  ,vhi.y+1,vlo.z  ,icomp)
                    + vel(vhi.x  ,vhi.y  ,vlo.z-1,icomp)
                    - vel(vhi.x  ,vhi.y  ,vlo.z  ,icomp) * 2.0;
            } else if (vhi.x == dhi.x and vhi.y == dhi.y) {
                vel      (vhi.x+1,vhi.y+1,vlo.z-1,icomp)
                    = vel(vhi.x+1,vhi.y  ,vlo.z-1,icomp)
                    + vel(vhi.x  ,vhi.y+1,vlo.z-1,icomp)
                    - vel(vhi.x  ,vhi.y  ,vlo.z-1,icomp);
            } else if (vhi.x == dhi.x and vlo.z == dlo.z) {
                vel      (vhi.x+1,vhi.y+1,vlo.z-1,icomp)
                    = vel(vhi.x+1,vhi.y+1,vlo.z  ,icomp)
                    + vel(vhi.x  ,vhi.y+1,vlo.z-1,icomp)
                    - vel(vhi.x  ,vhi.y+1,vlo.z  ,icomp);
            } else if (vhi.y == dhi.y and vlo.z == dlo.z) {
                vel      (vhi.x+1,vhi.y+1,vlo.z-1,icomp)
                    = vel(vhi.x+1,vhi.y+1,vlo.z  ,icomp)
                    + vel(vhi.x+1,vhi.y  ,vlo.z-1,icomp)
                    - vel(vhi.x+1,vhi.y  ,vlo.z  ,icomp);
            } else if (vhi.x == dhi.x) {
                int offset = AMREX_SPACEDIM * oxhi;
                mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                   vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
            } else if (vhi.y == dhi.y) {
                int offset = AMREX_SPACEDIM * oyhi;
                mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                   vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
            } else if (vlo.z == dlo.z) {
                int offset = AMREX_SPACEDIM * ozlo;
                mllinop_apply_bc_z(Orientation::low, bx, blen.z,
                                   vel, mzlo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalzlo, maxorder, dxinv[zdir], inhomog, icomp);
            } else if (mxhi(vhi.x+1,vhi.y+1,vlo.z-1) != BndryData::covered) {
                if (myhi(vhi.x,vhi.y+1,vlo.z-1) == BndryData::covered) {
                    int offset = AMREX_SPACEDIM * oxhi;
                    mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                       vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
                } else if (mxhi(vhi.x+1,vhi.y,vlo.z-1) == BndryData::covered) {
                    int offset = AMREX_SPACEDIM * oyhi;
                    mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                       vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
                } else {
                    int offset = AMREX_SPACEDIM * ozlo;
                    mllinop_apply_bc_z(Orientation::low, bx, blen.z,
                                       vel, mzlo, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalzlo, maxorder, dxinv[zdir], inhomog, icomp);
                }
            }
            break;
        }
        case 4: {
            // xlo & ylo & zhi
            Box bx = amrex::adjCellHi(amrex::adjCellLo(amrex::adjCellLo(vbox,xdir,1),ydir,1),zdir,1);
            if (vlo.x == dlo.x and vlo.y == dlo.y and vhi.z == dhi.z) {
                vel      (vlo.x-1, vlo.y-1, vhi.z+1,icomp)
                    = vel(vlo.x-1, vlo.y  , vhi.z  ,icomp)
                    + vel(vlo.x  , vlo.y-1, vhi.z  ,icomp)
                    + vel(vlo.x  , vlo.y  , vhi.z+1,icomp)
                    - vel(vlo.x  , vlo.y  , vhi.z  ,icomp) * 2.0;
            } else if (vlo.x == dlo.x and vlo.y == dlo.y) {
                vel      (vlo.x-1, vlo.y-1, vhi.z+1,icomp)
                    = vel(vlo.x-1, vlo.y  , vhi.z+1,icomp)
                    + vel(vlo.x  , vlo.y-1, vhi.z+1,icomp)
                    - vel(vlo.x  , vlo.y  , vhi.z+1,icomp);
            } else if (vlo.x == dlo.x and vhi.z == dhi.z) {
                vel      (vlo.x-1, vlo.y-1, vhi.z+1,icomp)
                    = vel(vlo.x-1, vlo.y-1, vhi.z  ,icomp)
                    + vel(vlo.x  , vlo.y-1, vhi.z+1,icomp)
                    - vel(vlo.x  , vlo.y-1, vhi.z  ,icomp);
            } else if (vlo.y == dlo.y and vhi.z == dhi.z) {
                vel      (vlo.x-1, vlo.y-1, vhi.z+1,icomp)
                    = vel(vlo.x-1, vlo.y-1, vhi.z  ,icomp)
                    + vel(vlo.x-1, vlo.y  , vhi.z+1,icomp)
                    - vel(vlo.x-1, vlo.y  , vhi.z  ,icomp);
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
            } else if (vhi.z == dhi.z) {
                int offset = AMREX_SPACEDIM * ozhi;
                mllinop_apply_bc_z(Orientation::high, bx, blen.z,
                                   vel, mzhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalzhi, maxorder, dxinv[zdir], inhomog, icomp);
            } else if (mxlo(vlo.x-1,vlo.y-1,vhi.z+1) != BndryData::covered) {
                if (mylo(vlo.x,vlo.y-1,vhi.z+1) == BndryData::covered) {
                    int offset = AMREX_SPACEDIM * oxlo;
                    mllinop_apply_bc_x(Orientation::low, bx, blen.x,
                                       vel, mxlo, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalxlo, maxorder, dxinv[xdir], inhomog, icomp);
                } else if (mxlo(vlo.x-1,vlo.y,vhi.z+1) == BndryData::covered) {
                    int offset = AMREX_SPACEDIM * oylo;
                    mllinop_apply_bc_y(Orientation::low, bx, blen.y,
                                       vel, mylo, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalylo, maxorder, dxinv[ydir], inhomog, icomp);
                } else {
                    int offset = AMREX_SPACEDIM * ozhi;
                    mllinop_apply_bc_z(Orientation::high, bx, blen.z,
                                       vel, mzhi, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalzhi, maxorder, dxinv[zdir], inhomog, icomp);
                }
            }
            break;
        }
        case 5: {
            // xhi & ylo & zhi
            Box bx = amrex::adjCellHi(amrex::adjCellLo(amrex::adjCellHi(vbox,xdir,1),ydir,1),zdir,1);
            if (vhi.x == dhi.x and vlo.y == dlo.y and vhi.z == dhi.z) {
                vel      (vhi.x+1,vlo.y-1,vhi.z+1,icomp)
                    = vel(vhi.x+1,vlo.y  ,vhi.z  ,icomp)
                    + vel(vhi.x  ,vlo.y-1,vhi.z  ,icomp)
                    + vel(vhi.x  ,vlo.y  ,vhi.z+1,icomp)
                    - vel(vhi.x  ,vlo.y  ,vhi.z  ,icomp) * 2.0;
            } else if (vhi.x == dhi.x and vlo.y == dlo.y) {
                vel      (vhi.x+1,vlo.y-1,vhi.z+1,icomp)
                    = vel(vhi.x+1,vlo.y  ,vhi.z+1,icomp)
                    + vel(vhi.x  ,vlo.y-1,vhi.z+1,icomp)
                    - vel(vhi.x  ,vlo.y  ,vhi.z+1,icomp);
            } else if (vhi.x == dhi.x and vhi.z == dhi.z) {
                vel      (vhi.x+1,vlo.y-1,vhi.z+1,icomp)
                    = vel(vhi.x+1,vlo.y-1,vhi.z  ,icomp)
                    + vel(vhi.x  ,vlo.y-1,vhi.z+1,icomp)
                    - vel(vhi.x  ,vlo.y-1,vhi.z  ,icomp);
            } else if (vlo.y == dlo.y and vhi.z == dhi.z) {
                vel      (vhi.x+1,vlo.y-1,vhi.z+1,icomp)
                    = vel(vhi.x+1,vlo.y-1,vhi.z  ,icomp)
                    + vel(vhi.x+1,vlo.y  ,vhi.z+1,icomp)
                    - vel(vhi.x+1,vlo.y  ,vhi.z  ,icomp);
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
            } else if (vhi.z == dhi.z) {
                int offset = AMREX_SPACEDIM * ozhi;
                mllinop_apply_bc_z(Orientation::high, bx, blen.z,
                                   vel, mzhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalzhi, maxorder, dxinv[zdir], inhomog, icomp);
            } else  if (mxhi(vhi.x+1,vlo.y-1,vhi.z+1) != BndryData::covered) {
                if (mylo(vhi.x,vlo.y-1,vhi.z+1) == BndryData::covered) {
                    int offset = AMREX_SPACEDIM * oxhi;
                    mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                       vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
                } else if (mxhi(vhi.x+1,vlo.y,vhi.z+1) == BndryData::covered) {
                    int offset = AMREX_SPACEDIM * oylo;
                    mllinop_apply_bc_y(Orientation::low, bx, blen.y,
                                       vel, mylo, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalylo, maxorder, dxinv[ydir], inhomog, icomp);
                } else {
                    int offset = AMREX_SPACEDIM * ozhi;
                    mllinop_apply_bc_z(Orientation::high, bx, blen.z,
                                       vel, mzhi, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalzhi, maxorder, dxinv[zdir], inhomog, icomp);
                }
            }
            break;
        }
        case 6: {
            // xlo & yhi & zhi
            Box bx = amrex::adjCellHi(amrex::adjCellHi(amrex::adjCellLo(vbox,xdir,1),ydir,1),zdir,1);
            if (vlo.x == dlo.x and vhi.y == dhi.y and vhi.z == dhi.z) {
                vel      (vlo.x-1,vhi.y+1,vhi.z+1,icomp)
                    = vel(vlo.x-1,vhi.y  ,vhi.z  ,icomp)
                    + vel(vlo.x  ,vhi.y+1,vhi.z  ,icomp)
                    + vel(vlo.x  ,vhi.y  ,vhi.z+1,icomp)
                    - vel(vlo.x  ,vhi.y  ,vhi.z  ,icomp) * 2.0;
            } else if (vlo.x == dlo.x and vhi.y == dhi.y) {
                vel      (vlo.x-1,vhi.y+1,vhi.z+1,icomp)
                    = vel(vlo.x-1,vhi.y  ,vhi.z+1,icomp)
                    + vel(vlo.x  ,vhi.y+1,vhi.z+1,icomp)
                    - vel(vlo.x  ,vhi.y  ,vhi.z+1,icomp);
            } else if (vlo.x == dlo.x and vhi.z == dhi.z) {
                vel      (vlo.x-1,vhi.y+1,vhi.z+1,icomp)
                    = vel(vlo.x-1,vhi.y+1,vhi.z  ,icomp)
                    + vel(vlo.x  ,vhi.y+1,vhi.z+1,icomp)
                    - vel(vlo.x  ,vhi.y+1,vhi.z  ,icomp);
            } else if (vhi.y == dhi.y and vhi.z == dhi.z) {
                vel      (vlo.x-1,vhi.y+1,vhi.z+1,icomp)
                    = vel(vlo.x-1,vhi.y+1,vhi.z  ,icomp)
                    + vel(vlo.x-1,vhi.y  ,vhi.z+1,icomp)
                    - vel(vlo.x-1,vhi.y  ,vhi.z  ,icomp);
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
            } else if (vhi.z == dhi.z) {
                int offset = AMREX_SPACEDIM * ozhi;
                mllinop_apply_bc_z(Orientation::high, bx, blen.z,
                                   vel, mzhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalzhi, maxorder, dxinv[zdir], inhomog, icomp);
            } else if (mxlo(vlo.x-1,vhi.y+1,vhi.z+1) != BndryData::covered) {
                if (myhi(vlo.x,vhi.y+1,vhi.z+1) == BndryData::covered) {
                    int offset = AMREX_SPACEDIM * oxlo;
                    mllinop_apply_bc_x(Orientation::low, bx, blen.x,
                                       vel, mxlo, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalxlo, maxorder, dxinv[xdir], inhomog, icomp);
                } else if (mxlo(vlo.x-1,vhi.y,vhi.z+1) == BndryData::covered) {
                    int offset = AMREX_SPACEDIM * oyhi;
                    mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                       vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
                } else {
                    int offset = AMREX_SPACEDIM * ozhi;
                    mllinop_apply_bc_z(Orientation::high, bx, blen.z,
                                       vel, mzhi, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalzhi, maxorder, dxinv[zdir], inhomog, icomp);
                }
            }
            break;
        }
        case 7: {
            // xhi & yhi & zhi
            Box bx = amrex::adjCellHi(amrex::adjCellHi(amrex::adjCellHi(vbox,xdir,1),ydir,1),zdir,1);
            if (vhi.x == dhi.x and vhi.y == dhi.y and vhi.z == dhi.z) {
                vel      (vhi.x+1,vhi.y+1,vhi.z+1,icomp)
                    = vel(vhi.x+1,vhi.y  ,vhi.z  ,icomp)
                    + vel(vhi.x  ,vhi.y+1,vhi.z  ,icomp)
                    + vel(vhi.x  ,vhi.y  ,vhi.z+1,icomp)
                    - vel(vhi.x  ,vhi.y  ,vhi.z  ,icomp) * 2.0;
            } else if (vhi.x == dhi.x and vhi.y == dhi.y) {
                vel      (vhi.x+1,vhi.y+1,vhi.z+1,icomp)
                    = vel(vhi.x+1,vhi.y  ,vhi.z+1,icomp)
                    + vel(vhi.x  ,vhi.y+1,vhi.z+1,icomp)
                    - vel(vhi.x  ,vhi.y  ,vhi.z+1,icomp);
            } else if (vhi.x == dhi.x and vhi.z == dhi.z) {
                vel      (vhi.x+1,vhi.y+1,vhi.z+1,icomp)
                    = vel(vhi.x+1,vhi.y+1,vhi.z  ,icomp)
                    + vel(vhi.x  ,vhi.y+1,vhi.z+1,icomp)
                    - vel(vhi.x  ,vhi.y+1,vhi.z  ,icomp);
            } else if (vhi.y == dhi.y and vhi.z == dhi.z) {
                vel      (vhi.x+1,vhi.y+1,vhi.z+1,icomp)
                    = vel(vhi.x+1,vhi.y+1,vhi.z  ,icomp)
                    + vel(vhi.x+1,vhi.y  ,vhi.z+1,icomp)
                    - vel(vhi.x+1,vhi.y  ,vhi.z  ,icomp);
            } else if (vhi.x == dhi.x) {
                int offset = AMREX_SPACEDIM * oxhi;
                mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                   vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
            } else if (vhi.y == dhi.y) {
                int offset = AMREX_SPACEDIM * oyhi;
                mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                   vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
            } else if (vhi.z == dhi.z) {
                int offset = AMREX_SPACEDIM * ozhi;
                mllinop_apply_bc_z(Orientation::high, bx, blen.z,
                                   vel, mzhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalzhi, maxorder, dxinv[zdir], inhomog, icomp);
            } else if (mxhi(vhi.x+1,vhi.y+1,vhi.z+1) != BndryData::covered) {
                if (myhi(vhi.x,vhi.y+1,vhi.z+1) == BndryData::covered) {
                    int offset = AMREX_SPACEDIM * oxhi;
                    mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                       vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
                } else if (mxhi(vhi.x+1,vhi.y,vhi.z+1) == BndryData::covered) {
                    int offset = AMREX_SPACEDIM * oyhi;
                    mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                       vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
                } else {
                    int offset = AMREX_SPACEDIM * ozhi;
                    mllinop_apply_bc_z(Orientation::high, bx, blen.z,
                                       vel, mzhi, bct[offset+icomp], bcl[offset+icomp],
                                       bcvalzhi, maxorder, dxinv[zdir], inhomog, icomp);
                }
            }
            break;
        }
        default: {}
        }
    }
}
  
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_edges (int iedge, Box const& vbox, // vbox: the valid box
                          Array4<Real> const& vel,
                          Array4<int const> const& mxlo,
                          Array4<int const> const& mylo,
                          Array4<int const> const& mzlo,
                          Array4<int const> const& mxhi,
                          Array4<int const> const& myhi,
                          Array4<int const> const& mzhi,
                          Array4<Real const> const& bcvalxlo,
                          Array4<Real const> const& bcvalylo,
                          Array4<Real const> const& bcvalzlo,
                          Array4<Real const> const& bcvalxhi,
                          Array4<Real const> const& bcvalyhi,
                          Array4<Real const> const& bcvalzhi,
                          GpuArray<BoundCond,2*AMREX_SPACEDIM*AMREX_SPACEDIM> const& bct,
                          GpuArray<Real,2*AMREX_SPACEDIM*AMREX_SPACEDIM> const& bcl,
			  int inhomog, int maxorder,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv, Box const& domain) noexcept
{
    constexpr int oxlo = 0;
    constexpr int oylo = 1;
    constexpr int ozlo = 2;
    constexpr int oxhi = 3;
    constexpr int oyhi = 4;
    constexpr int ozhi = 5;
    constexpr int xdir = 0;
    constexpr int ydir = 1;
    constexpr int zdir = 2;
    const auto blen = amrex::length(vbox);
    const auto vlo  = amrex::lbound(vbox);
    const auto vhi  = amrex::ubound(vbox);
    const auto dlo  = amrex::lbound(domain);
    const auto dhi  = amrex::ubound(domain);
    for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
        switch (iedge) {
        case 0: {
            // xlo & ylo
            if (vlo.x == dlo.x and vlo.y == dlo.y) {
                for (int k = vlo.z; k <= vhi.z; ++k) {
                    vel      (vlo.x-1,vlo.y-1,k,icomp)
                        = vel(vlo.x  ,vlo.y-1,k,icomp)
                        + vel(vlo.x-1,vlo.y  ,k,icomp)
                        - vel(vlo.x  ,vlo.y  ,k,icomp);
                }
            } else if (vlo.x == dlo.x) {
                Box bx = amrex::adjCellLo(amrex::adjCellLo(vbox,xdir,1),ydir,1);
                int offset = AMREX_SPACEDIM * oxlo;
                mllinop_apply_bc_x(Orientation::low, bx, blen.x,
                                   vel, mxlo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxlo, maxorder, dxinv[xdir], inhomog, icomp);
            } else if (vlo.y == dlo.y) {
                Box bx = amrex::adjCellLo(amrex::adjCellLo(vbox,xdir,1),ydir,1);
                int offset = AMREX_SPACEDIM * oylo;
                mllinop_apply_bc_y(Orientation::low, bx, blen.y,
                                   vel, mylo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalylo, maxorder, dxinv[ydir], inhomog, icomp);
            } else {
                for (int k = vlo.z; k <= vhi.z; ++k) {
                    if (mxlo(vlo.x-1,vlo.y-1,k) != BndryData::covered) {
                        Box bx(IntVect(vlo.x-1,vlo.y-1,k),IntVect(vlo.x-1,vlo.y-1,k));
                        if (mylo(vlo.x,vlo.y-1,k) == BndryData::covered) {
                            int offset = AMREX_SPACEDIM * oxlo;
                            mllinop_apply_bc_x(Orientation::low, bx, blen.x,
                                               vel, mxlo, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalxlo, maxorder, dxinv[xdir], inhomog, icomp);
                        } else {
                            int offset = AMREX_SPACEDIM * oylo;
                            mllinop_apply_bc_y(Orientation::low, bx, blen.y,
                                               vel, mylo, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalylo, maxorder, dxinv[ydir], inhomog, icomp);
                        }
                    }
                }
            }
            break;
        }
        case 1: {
            // xhi & ylo
            if (vhi.x == dhi.x and vlo.y == dlo.y) {
                for (int k = vlo.z; k <= vhi.z; ++k) {
                    vel      (vhi.x+1,vlo.y-1,k,icomp)
                        = vel(vhi.x  ,vlo.y-1,k,icomp)
                        + vel(vhi.x+1,vlo.y  ,k,icomp)
                        - vel(vhi.x  ,vlo.y  ,k,icomp);
                }
            } else if (vhi.x == dhi.x) {
                Box bx = amrex::adjCellLo(amrex::adjCellHi(vbox,xdir,1),ydir,1);
                int offset = AMREX_SPACEDIM * oxhi;
                mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                   vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
            } else if (vlo.y == dlo.y) {
                Box bx = amrex::adjCellLo(amrex::adjCellHi(vbox,xdir,1),ydir,1);
                int offset = AMREX_SPACEDIM * oylo;
                mllinop_apply_bc_y(Orientation::low, bx, blen.y,
                                   vel, mylo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalylo, maxorder, dxinv[ydir], inhomog, icomp);
            } else {
                for (int k = vlo.z; k <= vhi.z; ++k) {
                    if (mxhi(vhi.x+1,vlo.y-1,k) != BndryData::covered) {
                        Box bx(IntVect(vhi.x+1,vlo.y-1,k),IntVect(vhi.x+1,vlo.y-1,k));
                        if (mylo(vhi.x,vlo.y-1,k) == BndryData::covered) {
                            int offset = AMREX_SPACEDIM * oxhi;
                            mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                               vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
                        } else {
                            int offset = AMREX_SPACEDIM * oylo;
                            mllinop_apply_bc_y(Orientation::low, bx, blen.y,
                                               vel, mylo, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalylo, maxorder, dxinv[ydir], inhomog, icomp);
                        }
                    }
                }
            }
            break;
        }
        case 2: {
            // xlo & yhi
            if (vlo.x == dlo.x and vhi.y == dhi.y) {
                for (int k = vlo.z; k <= vhi.z; ++k) {
                    vel      (vlo.x-1,vhi.y+1,k,icomp)
                        = vel(vlo.x  ,vhi.y+1,k,icomp)
                        + vel(vlo.x-1,vhi.y  ,k,icomp)
                        - vel(vlo.x  ,vhi.y  ,k,icomp);
                }
            } else if (vlo.x == dlo.x) {
                Box bx = amrex::adjCellHi(amrex::adjCellLo(vbox,xdir,1),ydir,1);
                int offset = AMREX_SPACEDIM * oxlo;
                mllinop_apply_bc_x(Orientation::low, bx, blen.x,
                                   vel, mxlo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxlo, maxorder, dxinv[xdir], inhomog, icomp);
            } else if (vhi.y == dhi.y) {
                Box bx = amrex::adjCellHi(amrex::adjCellLo(vbox,xdir,1),ydir,1);
                int offset = AMREX_SPACEDIM * oyhi;
                mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                   vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
            } else {
                for (int k = vlo.z; k <= vhi.z; ++k) {
                    if (mxlo(vlo.x-1,vhi.y+1,k) != BndryData::covered) {
                        Box bx(IntVect(vlo.x-1,vhi.y+1,k),IntVect(vlo.x-1,vhi.y+1,k));
                        if (myhi(vlo.x,vhi.y+1,k) == BndryData::covered) {
                            int offset = AMREX_SPACEDIM * oxlo;
                            mllinop_apply_bc_x(Orientation::low, bx, blen.x,
                                               vel, mxlo, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalxlo, maxorder, dxinv[xdir], inhomog, icomp);
                        } else {
                            int offset = AMREX_SPACEDIM * oyhi;
                            mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                               vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
                        }
                    }
                }
            }
            break;
        }
        case 3: {
            // xhi & yhi
            if (vhi.x == dhi.x and vhi.y == dhi.y) {
                for (int k = vlo.z; k <= vhi.z; ++k) {
                    vel      (vhi.x+1,vhi.y+1,k,icomp)
                        = vel(vhi.x  ,vhi.y+1,k,icomp)
                        + vel(vhi.x+1,vhi.y  ,k,icomp)
                        - vel(vhi.x  ,vhi.y  ,k,icomp);
                }
            } else if (vhi.x == dhi.x) {
                Box bx = amrex::adjCellHi(amrex::adjCellHi(vbox,xdir,1),ydir,1);
                int offset = AMREX_SPACEDIM * oxhi;
                mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                   vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
            } else if (vhi.y == dhi.y) {
                Box bx = amrex::adjCellHi(amrex::adjCellHi(vbox,xdir,1),ydir,1);
                int offset = AMREX_SPACEDIM * oyhi;
                mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                   vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
            } else {
                for (int k = vlo.z; k <= vhi.z; ++k) {
                    if (mxhi(vhi.x+1,vhi.y+1,k) != BndryData::covered) {
                        Box bx(IntVect(vhi.x+1,vhi.y+1,k),IntVect(vhi.x+1,vhi.y+1,k));
                        if (myhi(vhi.x,vhi.y+1,k) == BndryData::covered) {
                            int offset = AMREX_SPACEDIM * oxhi;
                            mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                               vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
                        } else {
                            int offset = AMREX_SPACEDIM * oyhi;
                            mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                               vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
                        }
                    }
                }
            }
            break;
        }
        case 4: {
            // xlo & zlo
            if (vlo.x == dlo.x and vlo.z == dlo.z) {
                for (int j = vlo.y; j <= vhi.y; ++j) {
                    vel      (vlo.x-1,j,vlo.z-1,icomp)
                        = vel(vlo.x  ,j,vlo.z-1,icomp)
                        + vel(vlo.x-1,j,vlo.z  ,icomp)
                        - vel(vlo.x  ,j,vlo.z  ,icomp);
                }
            } else if (vlo.x == dlo.x) {
                Box bx = amrex::adjCellLo(amrex::adjCellLo(vbox,xdir,1),zdir,1);
                int offset = AMREX_SPACEDIM * oxlo;
                mllinop_apply_bc_x(Orientation::low, bx, blen.x,
                                   vel, mxlo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxlo, maxorder, dxinv[xdir], inhomog, icomp);
            } else if (vlo.z == dlo.z) {
                Box bx = amrex::adjCellLo(amrex::adjCellLo(vbox,xdir,1),zdir,1);
                int offset = AMREX_SPACEDIM * ozlo;
                mllinop_apply_bc_z(Orientation::low, bx, blen.z,
                                   vel, mzlo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalzlo, maxorder, dxinv[zdir], inhomog, icomp);
            } else {
                for (int j = vlo.y; j <= vhi.y; ++j) {
                    if (mxlo(vlo.x-1,j,vlo.z-1) != BndryData::covered) {
                        Box bx(IntVect(vlo.x-1,j,vlo.z-1),IntVect(vlo.x-1,j,vlo.z-1));
                        if (mzlo(vlo.x,j,vlo.z-1) == BndryData::covered) {
                            int offset = AMREX_SPACEDIM * oxlo;
                            mllinop_apply_bc_x(Orientation::low, bx, blen.x,
                                               vel, mxlo, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalxlo, maxorder, dxinv[xdir], inhomog, icomp);
                        } else {
                            int offset = AMREX_SPACEDIM * ozlo;
                            mllinop_apply_bc_z(Orientation::low, bx, blen.z,
                                               vel, mzlo, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalzlo, maxorder, dxinv[zdir], inhomog, icomp);
                        }
                    }
                }
            }
            break;
        }
        case 5: {
            // xhi & zlo
            if (vhi.x == dhi.x and vlo.z == dlo.z) {
                for (int j = vlo.y; j <= vhi.y; ++j) {
                    vel      (vhi.x+1,j,vlo.z-1,icomp)
                        = vel(vhi.x  ,j,vlo.z-1,icomp)
                        + vel(vhi.x+1,j,vlo.z  ,icomp)
                        - vel(vhi.x  ,j,vlo.z  ,icomp);
                }
            } else if (vhi.x == dhi.x) {
                Box bx = amrex::adjCellLo(amrex::adjCellHi(vbox,xdir,1),zdir,1);
                int offset = AMREX_SPACEDIM * oxhi;
                mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                   vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
            } else if (vlo.z == dlo.z) {
                Box bx = amrex::adjCellLo(amrex::adjCellHi(vbox,xdir,1),zdir,1);
                int offset = AMREX_SPACEDIM * ozlo;
                mllinop_apply_bc_z(Orientation::low, bx, blen.z,
                                   vel, mzlo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalzlo, maxorder, dxinv[zdir], inhomog, icomp);
            } else {
                for (int j = vlo.y; j <= vhi.y; ++j) {
                    if (mxhi(vhi.x+1,j,vlo.z-1) != BndryData::covered) {
                        Box bx(IntVect(vhi.x+1,j,vlo.z-1),IntVect(vhi.x+1,j,vlo.z-1));
                        if (mzlo(vhi.x,j,vlo.z-1) == BndryData::covered) {
                            int offset = AMREX_SPACEDIM * oxhi;
                            mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                               vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
                        } else {
                            int offset = AMREX_SPACEDIM * ozlo;
                            mllinop_apply_bc_z(Orientation::low, bx, blen.z,
                                               vel, mzlo, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalzlo, maxorder, dxinv[zdir], inhomog, icomp);
                        }
                    }
                }
            }
            break;
        }
        case 6: {
            // xlo & zhi
            if (vlo.x == dlo.x and vhi.z == dhi.z) {
                for (int j = vlo.y; j <= vhi.y; ++j) {
                    vel      (vlo.x-1,j,vhi.z+1,icomp)
                        = vel(vlo.x  ,j,vhi.z+1,icomp)
                        + vel(vlo.x-1,j,vhi.z  ,icomp)
                        - vel(vlo.x  ,j,vhi.z  ,icomp);
                }
            } else if (vlo.x == dlo.x) {
                Box bx = amrex::adjCellHi(amrex::adjCellLo(vbox,xdir,1),zdir,1);
                int offset = AMREX_SPACEDIM * oxlo;
                mllinop_apply_bc_x(Orientation::low, bx, blen.x,
                                   vel, mxlo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxlo, maxorder, dxinv[xdir], inhomog, icomp);
            } else if (vhi.z == dhi.z) {
                Box bx = amrex::adjCellHi(amrex::adjCellLo(vbox,xdir,1),zdir,1);
                int offset = AMREX_SPACEDIM * ozhi;
                mllinop_apply_bc_z(Orientation::high, bx, blen.z,
                                   vel, mzhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalzhi, maxorder, dxinv[zdir], inhomog, icomp);
            } else {
                for (int j = vlo.y; j <= vhi.y; ++j) {
                    if (mxlo(vlo.x-1,j,vhi.z+1) != BndryData::covered) {
                        Box bx(IntVect(vlo.x-1,j,vhi.z+1),IntVect(vlo.x-1,j,vhi.z+1));
                        if (mzhi(vlo.x,j,vhi.z+1) == BndryData::covered) {
                            int offset = AMREX_SPACEDIM * oxlo;
                            mllinop_apply_bc_x(Orientation::low, bx, blen.x,
                                               vel, mxlo, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalxlo, maxorder, dxinv[xdir], inhomog, icomp);
                        } else {
                            int offset = AMREX_SPACEDIM * ozhi;
                            mllinop_apply_bc_z(Orientation::high, bx, blen.z,
                                               vel, mzhi, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalzhi, maxorder, dxinv[zdir], inhomog, icomp);
                        }
                    }
                }
            }
            break;
        }
        case 7: {
            // xhi & zhi
            if (vhi.x == dhi.x and vhi.z == dhi.z) {
                for (int j = vlo.y; j <= vhi.y; ++j) {
                    vel      (vhi.x+1,j,vhi.z+1,icomp)
                        = vel(vhi.x  ,j,vhi.z+1,icomp)
                        + vel(vhi.x+1,j,vhi.z  ,icomp)
                        - vel(vhi.x  ,j,vhi.z  ,icomp);
                }
            } else if (vhi.x == dhi.x) {
                Box bx = amrex::adjCellHi(amrex::adjCellHi(vbox,xdir,1),zdir,1);
                int offset = AMREX_SPACEDIM * oxhi;
                mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                   vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
            } else if (vhi.z == dhi.z) {
                Box bx = amrex::adjCellHi(amrex::adjCellHi(vbox,xdir,1),zdir,1);
                int offset = AMREX_SPACEDIM * ozhi;
                mllinop_apply_bc_z(Orientation::high, bx, blen.z,
                                   vel, mzhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalzhi, maxorder, dxinv[zdir], inhomog, icomp);
            } else {
                for (int j = vlo.y; j <= vhi.y; ++j) {
                    if (mxhi(vhi.x+1,j,vhi.z+1) != BndryData::covered) {
                        Box bx(IntVect(vhi.x+1,j,vhi.z+1),IntVect(vhi.x+1,j,vhi.z+1));
                        if (mzhi(vhi.x,j,vhi.z+1) == BndryData::covered) {
                            int offset = AMREX_SPACEDIM * oxhi;
                            mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                               vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
                        } else {
                            int offset = AMREX_SPACEDIM * ozhi;
                            mllinop_apply_bc_z(Orientation::high, bx, blen.z,
                                               vel, mzhi, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalzhi, maxorder, dxinv[zdir], inhomog, icomp);
                        }
                    }
                }
            }
            break;
        }
        case 8: {
            // ylo & zlo
            if (vlo.y == dlo.y and vlo.z == dlo.z) {
                for (int i = vlo.x; i <= vhi.x; ++i) {
                    vel      (i,vlo.y-1,vlo.z-1,icomp)
                        = vel(i,vlo.y  ,vlo.z-1,icomp)
                        + vel(i,vlo.y-1,vlo.z  ,icomp)
                        - vel(i,vlo.y  ,vlo.z  ,icomp);
                }
            } else if (vlo.y == dlo.y) {
                Box bx = amrex::adjCellLo(amrex::adjCellLo(vbox,ydir,1),zdir,1);
                int offset = AMREX_SPACEDIM * oylo;
                mllinop_apply_bc_y(Orientation::low, bx, blen.y,
                                   vel, mylo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalylo, maxorder, dxinv[ydir], inhomog, icomp);
            } else if (vlo.z == dlo.z) {
                Box bx = amrex::adjCellLo(amrex::adjCellLo(vbox,ydir,1),zdir,1);
                int offset = AMREX_SPACEDIM * ozlo;
                mllinop_apply_bc_z(Orientation::low, bx, blen.z,
                                   vel, mzlo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalzlo, maxorder, dxinv[zdir], inhomog, icomp);
            } else {
                for (int i = vlo.x; i <= vhi.x; ++i) {
                    if (mylo(i,vlo.y-1,vlo.z-1) != BndryData::covered) {
                        Box bx(IntVect(i,vlo.y-1,vlo.z-1),IntVect(i,vlo.y-1,vlo.z-1));
                        if (mzlo(i,vlo.y,vlo.z-1) == BndryData::covered) {
                            int offset = AMREX_SPACEDIM * oylo;
                            mllinop_apply_bc_y(Orientation::low, bx, blen.y,
                                               vel, mylo, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalylo, maxorder, dxinv[ydir], inhomog, icomp);
                        } else {
                            int offset = AMREX_SPACEDIM * ozlo;
                            mllinop_apply_bc_z(Orientation::low, bx, blen.z,
                                               vel, mzlo, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalzlo, maxorder, dxinv[zdir], inhomog, icomp);
                        }
                    }
                }
            }
            break;
        }
        case 9: {
            // yhi & zlo
            if (vhi.y == dhi.y and vlo.z == dlo.z) {
                for (int i = vlo.x; i <= vhi.x; ++i) {
                    vel      (i,vhi.y+1,vlo.z-1,icomp)
                        = vel(i,vhi.y  ,vlo.z-1,icomp)
                        + vel(i,vhi.y+1,vlo.z  ,icomp)
                        - vel(i,vhi.y  ,vlo.z  ,icomp);
                }
            } else if (vhi.y == dhi.y) {
                Box bx = amrex::adjCellLo(amrex::adjCellHi(vbox,ydir,1),zdir,1);
                int offset = AMREX_SPACEDIM * oyhi;
                mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                   vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
            } else if (vlo.z == dlo.z) {
                Box bx = amrex::adjCellLo(amrex::adjCellHi(vbox,ydir,1),zdir,1);
                int offset = AMREX_SPACEDIM * ozlo;
                mllinop_apply_bc_z(Orientation::low, bx, blen.z,
                                   vel, mzlo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalzlo, maxorder, dxinv[zdir], inhomog, icomp);
            } else {
                for (int i = vlo.x; i <= vhi.x; ++i) {
                    if (myhi(i,vhi.y+1,vlo.z-1) != BndryData::covered) {
                        Box bx(IntVect(i,vhi.y+1,vlo.z-1),IntVect(i,vhi.y+1,vlo.z-1));
                        if (mzlo(i,vhi.y,vlo.z-1) == BndryData::covered) {
                            int offset = AMREX_SPACEDIM * oyhi;
                            mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                               vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
                        } else {
                            int offset = AMREX_SPACEDIM * ozlo;
                            mllinop_apply_bc_z(Orientation::low, bx, blen.z,
                                               vel, mzlo, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalzlo, maxorder, dxinv[zdir], inhomog, icomp);
                        }
                    }
                }
            }
            break;
        }
        case 10: {
            // ylo & zhi
            if (vlo.y == dlo.y and vhi.z == dhi.z) {
                for (int i = vlo.x; i <= vhi.x; ++i) {
                    vel      (i,vlo.y-1,vhi.z+1,icomp)
                        = vel(i,vlo.y  ,vhi.z+1,icomp)
                        + vel(i,vlo.y-1,vhi.z  ,icomp)
                        - vel(i,vlo.y  ,vhi.z  ,icomp);
                }
            } else if (vlo.y == dlo.y) {
                Box bx = amrex::adjCellHi(amrex::adjCellLo(vbox,ydir,1),zdir,1);
                int offset = AMREX_SPACEDIM * oylo;
                mllinop_apply_bc_y(Orientation::low, bx, blen.y,
                                   vel, mylo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalylo, maxorder, dxinv[ydir], inhomog, icomp);
            } else if (vhi.z == dhi.z) {
                Box bx = amrex::adjCellHi(amrex::adjCellLo(vbox,ydir,1),zdir,1);
                int offset = AMREX_SPACEDIM * ozhi;
                mllinop_apply_bc_z(Orientation::high, bx, blen.z,
                                   vel, mzhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalzhi, maxorder, dxinv[zdir], inhomog, icomp);
            } else {
                for (int i = vlo.x; i <= vhi.x; ++i) {
                    if (mylo(i,vlo.y-1,vhi.z+1) != BndryData::covered) {
                        Box bx(IntVect(i,vlo.y-1,vhi.z+1),IntVect(i,vlo.y-1,vhi.z+1));
                        if (mzhi(i,vlo.y,vhi.z+1) == BndryData::covered) {
                            int offset = AMREX_SPACEDIM * oylo;
                            mllinop_apply_bc_y(Orientation::low, bx, blen.y,
                                               vel, mylo, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalylo, maxorder, dxinv[ydir], inhomog, icomp);
                        } else {
                            int offset = AMREX_SPACEDIM * ozhi;
                            mllinop_apply_bc_z(Orientation::high, bx, blen.z,
                                               vel, mzhi, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalzhi, maxorder, dxinv[zdir], inhomog, icomp);
                        }
                    }
                }
            }
            break;
        }
        case 11: {
            // yhi & zhi
            if (vhi.y == dhi.y and vhi.z == dhi.z) {
                for (int i = vlo.x; i <= vhi.x; ++i) {
                    vel      (i,vhi.y+1,vhi.z+1,icomp)
                        = vel(i,vhi.y  ,vhi.z+1,icomp)
                        + vel(i,vhi.y+1,vhi.z  ,icomp)
                        - vel(i,vhi.y  ,vhi.z  ,icomp);
                }
            } else if (vhi.y == dhi.y) {
                Box bx = amrex::adjCellHi(amrex::adjCellHi(vbox,ydir,1),zdir,1);
                int offset = AMREX_SPACEDIM * oyhi;
                mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                   vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
            } else if (vhi.z == dhi.z) {
                Box bx = amrex::adjCellHi(amrex::adjCellHi(vbox,ydir,1),zdir,1);
                int offset = AMREX_SPACEDIM * ozhi;
                mllinop_apply_bc_z(Orientation::high, bx, blen.z,
                                   vel, mzhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalzhi, maxorder, dxinv[zdir], inhomog, icomp);
            } else {
                for (int i = vlo.x; i <= vhi.x; ++i) {
                    if (myhi(i,vhi.y+1,vhi.z+1) != BndryData::covered) {
                        Box bx(IntVect(i,vhi.y+1,vhi.z+1),IntVect(i,vhi.y+1,vhi.z+1));
                        if (mzhi(i,vhi.y,vhi.z+1) == BndryData::covered) {
                            int offset = AMREX_SPACEDIM * oyhi;
                            mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                               vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
                        } else {
                            int offset = AMREX_SPACEDIM * ozhi;
                            mllinop_apply_bc_z(Orientation::high, bx, blen.z,
                                               vel, mzhi, bct[offset+icomp], bcl[offset+icomp],
                                               bcvalzhi, maxorder, dxinv[zdir], inhomog, icomp);
                        }
                    }
                }
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
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dudy = (vel(i,j+1,k,0)+vel(i-1,j+1,k,0)-vel(i,j-1,k,0)-vel(i-1,j-1,k,0))*(0.25*dyi);
                Real dvdy = (vel(i,j+1,k,1)+vel(i-1,j+1,k,1)-vel(i,j-1,k,1)-vel(i-1,j-1,k,1))*(0.25*dyi);
                Real dudz = (vel(i,j,k+1,0)+vel(i-1,j,k+1,0)-vel(i,j,k-1,0)-vel(i-1,j,k-1,0))*(0.25*dzi);
                Real dwdz = (vel(i,j,k+1,2)+vel(i-1,j,k+1,2)-vel(i,j,k-1,2)-vel(i-1,j,k-1,2))*(0.25*dzi);
                Real divu = dvdy + dwdz;
                Real xif = kapx(i,j,k);
                Real mun = 0.75*(etax(i,j,k,0)-xif);  // restore the original eta
                Real mut =       etax(i,j,k,1);
                fx(i,j,k,0) = -mun*(-twoThirds*divu) - xif*divu;
                fx(i,j,k,1) = -mut*(dudy);
                fx(i,j,k,2) = -mut*(dudz);
            }
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
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dudx = (vel(i+1,j,k,0)+vel(i+1,j-1,k,0)-vel(i-1,j,k,0)-vel(i-1,j-1,k,0))*(0.25*dxi);
                Real dvdx = (vel(i+1,j,k,1)+vel(i+1,j-1,k,1)-vel(i-1,j,k,1)-vel(i-1,j-1,k,1))*(0.25*dxi);
                Real dvdz = (vel(i,j,k+1,1)+vel(i,j-1,k+1,1)-vel(i,j,k-1,1)-vel(i,j-1,k-1,1))*(0.25*dzi);
                Real dwdz = (vel(i,j,k+1,2)+vel(i,j-1,k+1,2)-vel(i,j,k-1,2)-vel(i,j-1,k-1,2))*(0.25*dzi);
                Real divu = dudx + dwdz;
                Real xif = kapy(i,j,k);
                Real mun = 0.75*(etay(i,j,k,1)-xif);  // restore the original eta
                Real mut =       etay(i,j,k,0);
                fy(i,j,k,0) = -mut*(dvdx);
                fy(i,j,k,1) = -mun*(-twoThirds*divu) - xif*divu;
                fy(i,j,k,2) = -mut*(dvdz);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms_fz (Box const& box, Array4<Real> const& fz,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& etaz,
                              Array4<Real const> const& kapz,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dudx = (vel(i+1,j,k,0)+vel(i+1,j,k-1,0)-vel(i-1,j,k,0)-vel(i-1,j,k-1,0))*(0.25*dxi);
                Real dwdx = (vel(i+1,j,k,2)+vel(i+1,j,k-1,2)-vel(i-1,j,k,2)-vel(i-1,j,k-1,2))*(0.25*dxi);
                Real dvdy = (vel(i,j+1,k,1)+vel(i,j+1,k-1,1)-vel(i,j-1,k,1)-vel(i,j-1,k-1,1))*(0.25*dyi);
                Real dwdy = (vel(i,j+1,k,2)+vel(i,j+1,k-1,2)-vel(i,j-1,k,2)-vel(i,j-1,k-1,2))*(0.25*dyi);
                Real divu = dudx + dvdy;
                Real xif = kapz(i,j,k);
                Real mun = 0.75*(etaz(i,j,k,2)-xif);  // restore the original eta
                Real mut =       etaz(i,j,k,0);
                fz(i,j,k,0) = -mut*(dwdx);
                fz(i,j,k,1) = -mut*(dwdy);
                fz(i,j,k,2) = -mun*(-twoThirds*divu) - xif*divu;
            }
        }
    } 
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms (Box const& box, Array4<Real> const& Ax,
                           Array4<Real const> const& fx,
                           Array4<Real const> const& fy,
                           Array4<Real const> const& fz,
                           GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                           Real bscalar) noexcept
{
    const Real dxi = bscalar * dxinv[0];
    const Real dyi = bscalar * dxinv[1];
    const Real dzi = bscalar * dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Ax(i,j,k,0) += dxi*(fx(i+1,j  ,k  ,0) - fx(i,j,k,0))
                    +          dyi*(fy(i  ,j+1,k  ,0) - fy(i,j,k,0))
                    +          dzi*(fz(i  ,j  ,k+1,0) - fz(i,j,k,0));
                Ax(i,j,k,1) += dxi*(fx(i+1,j  ,k  ,1) - fx(i,j,k,1))
                    +          dyi*(fy(i  ,j+1,k  ,1) - fy(i,j,k,1))
                    +          dzi*(fz(i  ,j  ,k+1,1) - fz(i,j,k,1));
                Ax(i,j,k,2) += dxi*(fx(i+1,j  ,k  ,2) - fx(i,j,k,2))
                    +          dyi*(fy(i  ,j+1,k  ,2) - fy(i,j,k,2))
                    +          dzi*(fz(i  ,j  ,k+1,2) - fz(i,j,k,2));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms_os (Box const& box, Array4<Real> const& Ax,
                              Array4<Real const> const& fx,
                              Array4<Real const> const& fy,
                              Array4<Real const> const& fz,
                              Array4<int const> const& osm,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Real bscalar) noexcept
{
    const Real dxi = bscalar * dxinv[0];
    const Real dyi = bscalar * dxinv[1];
    const Real dzi = bscalar * dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (osm(i,j,k) == 0) {
                    Ax(i,j,k,0) = 0.0;
                    Ax(i,j,k,1) = 0.0;
                    Ax(i,j,k,2) = 0.0;
                } else {
                    Ax(i,j,k,0) += dxi*(fx(i+1,j  ,k  ,0) - fx(i,j,k,0))
                        +          dyi*(fy(i  ,j+1,k  ,0) - fy(i,j,k,0))
                        +          dzi*(fz(i  ,j  ,k+1,0) - fz(i,j,k,0));
                    Ax(i,j,k,1) += dxi*(fx(i+1,j  ,k  ,1) - fx(i,j,k,1))
                        +          dyi*(fy(i  ,j+1,k  ,1) - fy(i,j,k,1))
                        +          dzi*(fz(i  ,j  ,k+1,1) - fz(i,j,k,1));
                    Ax(i,j,k,2) += dxi*(fx(i+1,j  ,k  ,2) - fx(i,j,k,2))
                        +          dyi*(fy(i  ,j+1,k  ,2) - fy(i,j,k,2))
                        +          dzi*(fz(i  ,j  ,k+1,2) - fz(i,j,k,2));
                }
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
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
              
                Real dudx = (vel(i,j,k,0) - vel(i-1,j,k,0))*dxi;
                Real dvdx = (vel(i,j,k,1) - vel(i-1,j,k,1))*dxi;
                Real dwdx = (vel(i,j,k,2) - vel(i-1,j,k,2))*dxi;
              
                Real dudy = (vel(i,j+1,k,0)+vel(i-1,j+1,k,0)-vel(i,j-1,k,0)-vel(i-1,j-1,k,0))*(0.25*dyi);
                Real dvdy = (vel(i,j+1,k,1)+vel(i-1,j+1,k,1)-vel(i,j-1,k,1)-vel(i-1,j-1,k,1))*(0.25*dyi);
                Real dwdy = (vel(i,j+1,k,2)+vel(i-1,j+1,k,2)-vel(i,j-1,k,2)-vel(i-1,j-1,k,2))*(0.25*dyi);
                
                Real dudz = (vel(i,j,k+1,0)+vel(i-1,j,k+1,0)-vel(i,j,k-1,0)-vel(i-1,j,k-1,0))*(0.25*dzi);
                Real dvdz = (vel(i,j,k+1,1)+vel(i-1,j,k+1,1)-vel(i,j,k-1,1)-vel(i-1,j,k-1,1))*(0.25*dzi);
                Real dwdz = (vel(i,j,k+1,2)+vel(i-1,j,k+1,2)-vel(i,j,k-1,2)-vel(i-1,j,k-1,2))*(0.25*dzi);
                
                fx(i,j,k,0) = dudx;
                fx(i,j,k,1) = dvdx;
                fx(i,j,k,2) = dwdx;
                fx(i,j,k,3) = dudy;
                fx(i,j,k,4) = dvdy;
                fx(i,j,k,5) = dwdy;
                fx(i,j,k,6) = dudz;
                fx(i,j,k,7) = dvdz;
                fx(i,j,k,8) = dwdz;
                
            }
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
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
              
                Real dudx = (vel(i+1,j,k,0)+vel(i+1,j-1,k,0)-vel(i-1,j,k,0)-vel(i-1,j-1,k,0))*(0.25*dxi);
                Real dvdx = (vel(i+1,j,k,1)+vel(i+1,j-1,k,1)-vel(i-1,j,k,1)-vel(i-1,j-1,k,1))*(0.25*dxi);
                Real dwdx = (vel(i+1,j,k,2)+vel(i+1,j-1,k,2)-vel(i-1,j,k,2)-vel(i-1,j-1,k,2))*(0.25*dxi);
                
                Real dudy = (vel(i,j,k,0) - vel(i,j-1,k,0))*dyi;
                Real dvdy = (vel(i,j,k,1) - vel(i,j-1,k,1))*dyi;
                Real dwdy = (vel(i,j,k,2) - vel(i,j-1,k,2))*dyi;
                
                Real dudz = (vel(i,j,k+1,0)+vel(i,j-1,k+1,0)-vel(i,j,k-1,0)-vel(i,j-1,k-1,0))*(0.25*dzi);
                Real dvdz = (vel(i,j,k+1,1)+vel(i,j-1,k+1,1)-vel(i,j,k-1,1)-vel(i,j-1,k-1,1))*(0.25*dzi);
                Real dwdz = (vel(i,j,k+1,2)+vel(i,j-1,k+1,2)-vel(i,j,k-1,2)-vel(i,j-1,k-1,2))*(0.25*dzi);
                
                fy(i,j,k,0) = dudx;
                fy(i,j,k,1) = dvdx;
                fy(i,j,k,2) = dwdx;
                fy(i,j,k,3) = dudy;
                fy(i,j,k,4) = dvdy;
                fy(i,j,k,5) = dwdy;
                fy(i,j,k,6) = dudz;
                fy(i,j,k,7) = dvdz;
                fy(i,j,k,8) = dwdz;

            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_vel_grads_fz (Box const& box, Array4<Real> const& fz,
                              Array4<Real const> const& vel,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
              
                Real dudx = (vel(i+1,j,k,0)+vel(i+1,j,k-1,0)-vel(i-1,j,k,0)-vel(i-1,j,k-1,0))*(0.25*dxi);
                Real dvdx = (vel(i+1,j,k,1)+vel(i+1,j,k-1,1)-vel(i-1,j,k,1)-vel(i-1,j,k-1,1))*(0.25*dxi);
                Real dwdx = (vel(i+1,j,k,2)+vel(i+1,j,k-1,2)-vel(i-1,j,k,2)-vel(i-1,j,k-1,2))*(0.25*dxi);
                
                Real dudy = (vel(i,j+1,k,0)+vel(i,j+1,k-1,0)-vel(i,j-1,k,0)-vel(i,j-1,k-1,0))*(0.25*dyi);
                Real dvdy = (vel(i,j+1,k,1)+vel(i,j+1,k-1,1)-vel(i,j-1,k,1)-vel(i,j-1,k-1,1))*(0.25*dyi);
                Real dwdy = (vel(i,j+1,k,2)+vel(i,j+1,k-1,2)-vel(i,j-1,k,2)-vel(i,j-1,k-1,2))*(0.25*dyi);
                
                Real dudz = (vel(i,j,k,0) - vel(i,j,k-1,0))*dzi;
                Real dvdz = (vel(i,j,k,1) - vel(i,j,k-1,1))*dzi;
                Real dwdz = (vel(i,j,k,2) - vel(i,j,k-1,2))*dzi;
                
                fz(i,j,k,0) = dudx;
                fz(i,j,k,1) = dvdx;
                fz(i,j,k,2) = dwdx;
                fz(i,j,k,3) = dudy;
                fz(i,j,k,4) = dvdy;
                fz(i,j,k,5) = dwdy;
                fz(i,j,k,6) = dudz;
                fz(i,j,k,7) = dvdz;
                fz(i,j,k,8) = dwdz;

            }
        }
    } 
}

}

#endif
