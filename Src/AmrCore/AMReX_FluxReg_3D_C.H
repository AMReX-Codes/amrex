#ifndef AMREX_FLUXREG_3D_C_H_
#define AMREX_FLUXREG_3D_C_H_

#include <AMReX_FArrayBox.H>

namespace amrex {

//
// Add fine grid flux to flux register.  Flux array is a fine grid
// edge based object, Register is a coarse grid edge based object.      
// It is assumed that the coarsened flux region contains the register
// region.
//
AMREX_GPU_HOST_DEVICE inline void
fluxreg_fineadd (Box const& bx, FArrayBox& regfab, const int rcomp,
                 FArrayBox const& flxfab, const int fcomp, const int ncomp,
                 const int dir, IntVect const& ratio, const Real mult)
{
    const auto lo  = amrex::lbound(bx);
    const auto hi  = amrex::ubound(bx);
    const auto reg = regfab.view(lo);

    const auto flo = amrex::lbound(flxfab.box());
    const auto flx = flxfab.view(flo);

    switch (dir) {
    case 0:
    {
        const int ic = lo.x;
        const int i = ic*ratio[0];
        for (int n = 0; n < ncomp; ++n) {
            for (int kc = lo.z; kc <= hi.z; ++kc) {
                for (int koff = 0; koff < ratio[2]; ++koff) {
                    const int k = ratio[2]*kc + koff;
                    for (int jc = lo.y; jc <= hi.y; ++jc) {
                        for (int joff = 0; joff < ratio[1]; ++joff) {
                            const int j = ratio[1]*jc + joff;
                            reg(ic-lo.x,jc-lo.y,kc-lo.z,n) += mult*flx(i-flo.x,j-flo.y,k-flo.z,n);
                        }
                    }
                }
            }
        }
        break;
    }
    case 1:
    {
        const int jc = lo.y;
        const int j = jc*ratio[1];
        for (int n = 0; n < ncomp; ++n) {
            for (int kc = lo.z; kc <= hi.z; ++kc) {
                for (int koff = 0; koff < ratio[2]; ++koff) {
                    const int k = ratio[2]*kc + koff;
                    for (int ic = lo.x; ic <= hi.x; ++ic) {
                        for (int ioff = 0; ioff < ratio[0]; ++ioff) {
                            const int i = ratio[0]*ic + ioff;
                            reg(ic-lo.x,jc-lo.y,kc-lo.z,n) += mult*flx(i-flo.x,j-flo.y,k-flo.z,n);
                        }
                    }
                }
            }
        }
        break;
    }
    default:
    {
        const int kc = lo.z;
        const int k = kc*ratio[2];
        for (int n = 0; n < ncomp; ++n) {
            for (int jc = lo.y; jc <= hi.y; ++jc) {
                for (int joff = 0; joff < ratio[1]; ++joff) {
                    const int j = ratio[1]*jc + joff;
                    for (int ic = lo.x; ic <= hi.x; ++ic) {
                        for (int ioff = 0; ioff < ratio[0]; ++ioff) {
                            const int i = ratio[0]*ic + ioff;
                            reg(ic-lo.x,jc-lo.y,kc-lo.z,n) += mult*flx(i-flo.x,j-flo.y,k-flo.z,n);
                        }
                    }
                }
            }
        }
    }
    }
}

//
// Add fine grid flux times area to flux register.  Flux array is a fine grid
// edge based object, Register is a coarse grid edge based object.      
// It is assumed that the coarsened flux region contains the register
// region.
//

}

#endif
