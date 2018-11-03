
#include <AMReX_FArrayBox.H>

using namespace amrex;

void flux_to_dudt_c (Box const& bx,
                     FArrayBox& dudtfab,
                     FArrayBox const& fxfab,
                     FArrayBox const& fyfab,
                     FArrayBox const& fzfab,
                     Array<Real,AMREX_SPACEDIM> const& dxinv)
{
    const auto len = bx.length3d();
    const int ncomp = dudtfab.nComp();
    const auto dudt = dudtfab.view(bx);
    const auto fx   =   fxfab.view(bx);
    const auto fy   =   fyfab.view(bx);
    const auto fz   =   fzfab.view(bx);

    for (int n = 0; n < 5; ++n) {
        for         (int k = 0; k < len[2]; ++k) {
            for     (int j = 0; j < len[1]; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len[0]; ++i) {
                    dudt(i,j,k,n) = dxinv[0] * (fx(i,j,k,n) - fx(i+1,j,k,n))
                        +           dxinv[1] * (fy(i,j,k,n) - fy(i,j+1,k,n))
                        +           dxinv[2] * (fz(i,j,k,n) - fz(i,j,k+1,n));
                }
            }
        }
    }
}
