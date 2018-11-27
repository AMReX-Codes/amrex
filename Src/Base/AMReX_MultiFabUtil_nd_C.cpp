
#include <AMReX_MultiFabUtil_C.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE
void amrex_int_to_real (Box const& bx, FArrayBox& rfab, IArrayBox const& ifab,
                        int rcomp, int icomp, int ncomp)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto rdata = rfab.view(lo,rcomp);
    const auto idata = ifab.view(lo,icomp);

    for (int n = 0; n < ncomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    rdata(i,j,k,n) = static_cast<Real>(idata(i,j,k,n));
                }
            }
        }
    }
}

//AMREX_GPU_HOST_DEVICE
//void amrex_fill_slice

}
