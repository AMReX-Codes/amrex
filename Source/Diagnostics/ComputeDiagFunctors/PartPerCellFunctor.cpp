#include "PartPerCellFunctor.H"
#include "WarpX.H"

using namespace amrex;

PartPerCellFunctor::PartPerCellFunctor(const amrex::MultiFab* mf_src, const int lev, const int ncomp)
    : ComputeDiagFunctor(ncomp), m_lev(lev)
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
    // Make alias MultiFab* pointing to the component of mf_dst where the
    // number of particles per cell is to be written.
    MultiFab * const mf_dst_dcomp = new MultiFab(mf_dst, amrex::make_alias, dcomp, 1);
    // Set value to 0, and increment the value in each cell with ppc.
    mf_dst_dcomp->setVal(0._rt);
    warpx.GetPartContainer().Increment(*mf_dst_dcomp, m_lev);
}
