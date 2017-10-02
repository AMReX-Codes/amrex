
#include <AmrCoreAdvPhysBC.H>
#include <AMReX_filcc_f.H>

AmrCoreAdvPhysBC::AmrCoreAdvPhysBC (const Geometry& geom, const BCRec& bcr)
    : m_geom(geom), m_bcr(bcr)
{ }

void
AmrCoreAdvPhysBC::define (const Geometry& geom, const BCRec& bcr)
{
    m_geom = geom;
    m_bcr = bcr;
}

void
AmrCoreAdvPhysBC::FillBoundary (MultiFab& mf, int dcomp, int ncomp, Real time)
{
    BL_PROFILE("AmrCoreAdvPhysBC::FillBoundary");

    if (mf.nGrow() == 0) return;
    
    if (m_geom.isAllPeriodic()) return;

    const Box&     domain      = m_geom.Domain();
    const Real*    dx          = m_geom.CellSize();
    const RealBox& prob_domain = m_geom.ProbDomain();
    const Real*    problo      = prob_domain.lo();

    // create a grown domain box containing valid + periodic cells
    Box gdomain = amrex::convert(domain, mf.boxArray().ixType());
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	if (m_geom.isPeriodic(i)) {
	    gdomain.grow(i, mf.nGrow());
	}
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = mf[mfi];
        const Box& fab_box = fab.box(); // including ghost cells
            
        // if the valid+ghost fab box contains cells not in the 
        // grown domain box, we need to fill physical ghost cells
        if (!gdomain.contains(fab_box))
        {
            amrex_fab_filcc(BL_TO_FORTRAN_FAB(fab),
                            BL_TO_FORTRAN_BOX(domain),
                            dx, problo,
                            m_bcr.data());
        }
    }

}


