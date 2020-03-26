#include "PartPerGridFunctor.H"
#include "Utils/Average.H"

using namespace amrex;

PartPerGridFunctor::PartPerGridFunctor(const amrex::MultiFab * const mf_src, const int lev, const int ncomp)
    : ComputeDiagFunctor(ncomp), m_lev(lev)
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
    // MultiFab containing number of particles per grid
    // (stored as constant for all cells in each grid)
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf_dst); mfi.isValid(); ++mfi) {
        mf_dst[mfi].setVal<RunOn::Host>(static_cast<Real>(npart_in_grid[mfi.index()]));
    }
}
