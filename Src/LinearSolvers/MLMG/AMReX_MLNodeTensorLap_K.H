#ifndef AMREX_MLNODETENSORLAP_K_H_
#define AMREX_MLNODETENSORLAP_K_H_
#include <AMReX_Config.H>

#include <AMReX_MLNodeTensorLaplacian.H>
#include <AMReX_LO_BCTYPES.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLNodeTensorLap_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLNodeTensorLap_2D_K.H>
#else
#include <AMReX_MLNodeTensorLap_3D_K.H>
#endif

namespace amrex {

// One color of a multi-color sweep: the nodes of bx whose index parities
// match offset.
template <unsigned int N>
void mlndtslap_gauss_seidel (Box const& bx, IntVect const& offset, Array4<Real> const& sol,
                             Array4<Real const> const& rhs, Array4<int const> const& msk,
                             GpuArray<Real,N> const& s) noexcept
{
    const auto lo = amrex::lbound(detail::strided_box(bx, IntVect(2), offset));
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; k += 2) {
    for (int j = lo.y; j <= hi.y; j += 2) {
    for (int i = lo.x; i <= hi.x; i += 2) {
        mlndtslap_gauss_seidel(i, j, k, sol, rhs, msk, s);
    }}}
}

}

#endif
