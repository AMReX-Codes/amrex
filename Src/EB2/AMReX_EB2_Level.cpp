
#include <AMReX_EB2_Level.H>

namespace amrex { namespace EB2 {

void
Level::fillVolFrac (MultiFab& vfrac, const Geometry& geom) const
{
    vfrac.setVal(1.0);
    vfrac.ParallelCopy(m_volfrac,0,0,1,0,vfrac.nGrow(),geom.periodicity());

    const std::vector<IntVect>& pshifts = geom.periodicity().shiftIntVect();

    Real cov_val = 0.0; // for covered cells
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector<std::pair<int,Box> > isects;
        for (MFIter mfi(vfrac); mfi.isValid(); ++mfi)
        {
            auto& fab = vfrac[mfi];
            const Box& bx = fab.box();
            for (const auto& iv : pshifts)
            {
                m_covered_grids.intersections(bx+iv, isects);
                for (const auto& is : isects) {
                    fab.setVal(cov_val, is.second-iv, 0, 1);
                }
            }
        }
    }
}

}}
