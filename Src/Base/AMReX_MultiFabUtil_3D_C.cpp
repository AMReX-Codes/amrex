
#include <AMReX_MultiFabUtil_C.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE
void amrex_avg_nd_to_cc (Box const& bx, FArrayBox& ccfab, FArrayBox const& ndfab,
                         int cccomp, int ndcomp, int ncomp)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto cc = ccfab.view(lo,cccomp);
    const auto nd = ndfab.view(lo,ndcomp);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = 0; k < len.z; ++k) {
        for (int j = 0; j < len.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            cc(i,j,k,n) = 0.125*( nd(i,j  ,k  ,n) + nd(i+1,j  ,k  ,n)
                                + nd(i,j+1,k  ,n) + nd(i+1,j+1,k  ,n)
                                + nd(i,j  ,k+1,n) + nd(i+1,j  ,k+1,n)
                                + nd(i,j+1,k+1,n) + nd(i+1,j+1,k+1,n));
        }}}
    }
}

}
