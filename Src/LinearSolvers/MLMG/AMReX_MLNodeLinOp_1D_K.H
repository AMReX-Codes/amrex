#ifndef AMREX_ML_NODE_LINOP_1D_K_H_
#define AMREX_ML_NODE_LINOP_1D_K_H_
#include <AMReX_Config.H>

namespace amrex {

template <typename T>
void mlndlap_bc_doit (Box const& vbx, Array4<T> const& a, Box const& domain,
                      GpuArray<bool,AMREX_SPACEDIM> const& bflo,
                      GpuArray<bool,AMREX_SPACEDIM> const& bfhi) noexcept
{
    Box gdomain = domain;
    int const idim = 0;
    if (! bflo[idim]) { gdomain.growLo(idim,1); }
    if (! bfhi[idim]) { gdomain.growHi(idim,1); }

    if (gdomain.strictly_contains(vbx)) { return; }

    const int offset = domain.cellCentered() ? 0 : 1;

    const auto dlo = domain.smallEnd(0);
    const auto dhi = domain.bigEnd(0);

    Box const& sbox = amrex::grow(vbx,1);
    AMREX_HOST_DEVICE_FOR_3D(sbox, i, j, k,
    {
        if (! gdomain.contains(IntVect(i))) {
            if (i == dlo-1 && bflo[0])
            {
                a(i,0,0) = a(i+1+offset, j, k);
            }
            else if (i == dhi+1 && bfhi[0])
            {
                a(i,0,0) = a(i-1-offset, j, k);
            }
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction (int i, int, int, Array4<Real> const& crse,
                          Array4<Real const> const& fine, Array4<int const> const& msk) noexcept
{
    int ii = i*2;
    if (msk(ii,0,0)) {
        crse(i,0,0) = Real(0.0);
    } else {
        crse(i,0,0) = Real(1./4.)  *(fine(ii-1,0,0)
                         + Real(2.)* fine(ii  ,0,0)
                         +           fine(ii+1,0,0));
    }
}

template <int rr>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction (int i, int, int, Array4<Real> const& crse,
                          Array4<Real const> const& fine, Array4<int const> const& msk,
                          Box const& fdom,
                          GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                          GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept

{
    const int ii = i*rr;
    if (msk(ii,0,0)) {
        crse(i,0,0) = Real(0.0);
    } else {
        const auto ndlo = fdom.smallEnd(0);
        const auto ndhi = fdom.bigEnd(0);
        Real tmp = Real(0.0);
        for (int ioff = -rr+1; ioff <= rr-1; ++ioff) {
            Real wx = rr - std::abs(ioff);
            int itmp = ii + ioff;
            if ((itmp < ndlo && (bclo[0] == LinOpBCType::Neumann ||
                                 bclo[0] == LinOpBCType::inflow)) ||
                (itmp > ndhi && (bchi[0] == LinOpBCType::Neumann ||
                                 bchi[0] == LinOpBCType::inflow))) {
                itmp = ii - ioff;
            }
            tmp += wx*fine(itmp,0,0);
        }
        crse(i,0,0) = tmp*(Real(1.0)/Real(rr*rr));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_restriction (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const&,
                          Array4<Real const> const&, Array4<int const> const&, int) noexcept
{
    amrex::Abort("mlndlap_semi_restriction: not implemented in 1D");
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_nodal_mask (int i, int, int, Array4<int> const& nmsk,
                             Array4<int const> const& cmsk) noexcept
{
    using namespace nodelap_detail;

    int s = cmsk(i-1,0,0) + cmsk(i,0,0);
    if (s == 2*crse_cell) {
        nmsk(i,0,0) = crse_node;
    } else if (s == 2*fine_cell) {
        nmsk(i,0,0) = fine_node;
    } else {
        nmsk(i,0,0) = crse_fine_node;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_dirichlet_mask (Box const& bx, Array4<int> const& dmsk,
                                 Array4<int const> const& omsk, Box const& dom,
                                 GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                                 GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    const auto lo = bx.smallEnd(0);
    const auto hi = bx.bigEnd(0);
    AMREX_PRAGMA_SIMD
    for (int i = lo; i <= hi; ++i) {
        if (!dmsk(i,0,0)) {
            dmsk(i,0,0) = (omsk(i-1,0,0) == 1 || omsk(i,0,0) == 1);
        }
    }

    const auto domlo = dom.smallEnd(0);
    const auto domhi = dom.bigEnd(0);

    if (bclo[0] == LinOpBCType::Dirichlet && lo == domlo) {
        dmsk(lo,0,0) = 1;
    }

    if (bchi[0] == LinOpBCType::Dirichlet && hi == domhi) {
        dmsk(hi,0,0) = 1;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_dot_mask (Box const& bx, Array4<Real> const& dmsk,
                           Array4<int const> const& omsk, Box const& dom,
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    const auto lo = bx.smallEnd(0);
    const auto hi = bx.bigEnd(0);

    AMREX_PRAGMA_SIMD
    for (int i = lo; i <= hi; ++i) {
        dmsk(i,0,0) = static_cast<Real>(omsk(i,0,0));
    }

    const auto domlo = dom.smallEnd(0);
    const auto domhi = dom.bigEnd(0);

    if ((bclo[0] == LinOpBCType::Neumann || bclo[0] == LinOpBCType::inflow)
        && lo == domlo)
    {
        dmsk(lo,0,0) *= Real(0.5);
    }

    if ((bchi[0] == LinOpBCType::Neumann || bchi[0] == LinOpBCType::inflow)
        && hi == domhi)
    {
        dmsk(hi,0,0) *= Real(0.5);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_dot_mask (Box const& bx, Array4<Real> const& dmsk,
                           Array4<int const> const& omsk,
                           Array4<int const> const& fmsk, Box const& dom,
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    const auto lo = bx.smallEnd(0);
    const auto hi = bx.bigEnd(0);

    AMREX_PRAGMA_SIMD
    for (int i = lo; i <= hi; ++i) {
        if (fmsk(i,0,0) == 0) {
            dmsk(i,0,0) = static_cast<Real>(omsk(i,0,0));
        } else {
            dmsk(i,0,0) = Real(0);
        }
    }

    const auto domlo = dom.smallEnd(0);
    const auto domhi = dom.bigEnd(0);

    if ((bclo[0] == LinOpBCType::Neumann || bclo[0] == LinOpBCType::inflow)
        && lo == domlo)
    {
        dmsk(lo,0,0) *= Real(0.5);
    }

    if ((bchi[0] == LinOpBCType::Neumann || bchi[0] == LinOpBCType::inflow)
        && hi == domhi)
    {
        dmsk(hi,0,0) *= Real(0.5);
    }
}

}

#endif
