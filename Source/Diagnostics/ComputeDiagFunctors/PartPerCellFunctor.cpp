#include "PartPerCellFunctor.H"
#include "WarpX.H"
#include "Utils/Average.H"

using namespace amrex;

PartPerCellFunctor::PartPerCellFunctor(const amrex::MultiFab* mf_src, const int lev, amrex::IntVect crse_ratio, const int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_lev(lev)
{
    // mf_src will not be used, let's make sure it's null.
    AMREX_ALWAYS_ASSERT(mf_src == nullptr);
    // Write only in one output component.
    AMREX_ALWAYS_ASSERT(ncomp == 1);
}

void
PartPerCellFunctor::operator()(amrex::MultiFab& mf_dst, const int dcomp) const
{
    auto& warpx = WarpX::GetInstance();
    // Guard cell is set to 1 for generality. However, for a cell-centered
    // output Multifab, mf_avg, the guard-cell data is not needed especially considering
    // the operations performend in the CoarsenAndInterpolate function.
    constexpr int ng = 1;
    // Temporary cell-centered, single-component MultiFab for storing particles per cell.
    MultiFab ppc_mf(warpx.boxArray(m_lev), warpx.DistributionMap(m_lev), 1, ng);
    // Set value to 0, and increment the value in each cell with ppc.
    ppc_mf.setVal(0._rt);
    // Compute ppc which includes a summation over all species.
    warpx.GetPartContainer().Increment(ppc_mf, m_lev);
    // Coarsen and interpolate from ppc_mf to the output diagnostic MultiFab, mf_dst.
    Average::CoarsenAndInterpolate(mf_dst, ppc_mf, dcomp, 0, nComp(), 0, m_crse_ratio);
}
