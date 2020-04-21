#include "PartPerGridFunctor.H"
#include "Utils/Average.H"

using namespace amrex;

PartPerGridFunctor::PartPerGridFunctor(const amrex::MultiFab * const mf_src, const int lev, const amrex::IntVect crse_ratio, const int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_lev(lev)
{
    // mf_src will not be used, let's make sure it's null.
    AMREX_ALWAYS_ASSERT(mf_src == nullptr);
    // Write only in one output component.
    AMREX_ALWAYS_ASSERT(ncomp == 1);
}

void
PartPerGridFunctor::operator()(amrex::MultiFab& mf_dst, const int dcomp) const
{
    auto& warpx = WarpX::GetInstance();
    const Vector<long>& npart_in_grid = warpx.GetPartContainer().NumberOfParticlesInGrid(m_lev);
    // Guard cell is set to 1 for generality. However, for a cell-centered
    // output Multifab, mf_dst, the guard-cell data is not needed especially considering
    // the operations performend in the CoarsenAndInterpolate function.
    constexpr int ng = 1;
    // Temporary MultiFab containing number of particles per grid.
    // (stored as constant for all cells in each grid)
    MultiFab ppg_mf(warpx.boxArray(m_lev), warpx.DistributionMap(m_lev), 1, ng);
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(ppg_mf); mfi.isValid(); ++mfi) {
        ppg_mf[mfi].setVal<RunOn::Host>(static_cast<Real>(npart_in_grid[mfi.index()]));
    }

    // Coarsen and interpolate from ppg_mf to the output diagnostic MultiFab, mf_dst.
    Average::CoarsenAndInterpolate(mf_dst, ppg_mf, dcomp, 0, nComp(), 0, m_crse_ratio);
}
