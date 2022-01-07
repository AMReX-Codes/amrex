#include "CNS_derive.H"
#include "CNS.H"
#include "CNS_parm.H"

using namespace amrex;

void cns_derpres (const Box& bx, FArrayBox& pfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& rhoefab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    auto const rhoe = rhoefab.array();
    auto       p    = pfab.array();
    Parm const* parm = CNS::d_parm;
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        p(i,j,k,dcomp) = (parm->eos_gamma-1.)*rhoe(i,j,k);
    });
}

void cns_dervel (const Box& bx, FArrayBox& velfab, int dcomp, int ncomp,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    amrex::ignore_unused(ncomp);
    AMREX_ASSERT(ncomp == AMREX_SPACEDIM);
    auto const dat = datfab.array();
    auto       vel = velfab.array();
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        AMREX_D_TERM(vel(i,j,k,dcomp  ) = dat(i,j,k,1)/dat(i,j,k,0);,
                     vel(i,j,k,dcomp+1) = dat(i,j,k,2)/dat(i,j,k,0);,
                     vel(i,j,k,dcomp+2) = dat(i,j,k,3)/dat(i,j,k,0);)
    });
}
