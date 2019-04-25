
#include <AMReX_AuxBoundaryData.H>
#include <AMReX_Print.H>

#ifdef BL_LAZY
#include <AMReX_Lazy.H>
#endif

namespace amrex {

AuxBoundaryData::AuxBoundaryData () noexcept
    :
    m_ngrow(0),
    m_empty(false),
    m_initialized(false)
{}

AuxBoundaryData::AuxBoundaryData (const BoxArray& ba,
                                  int             n_grow,
                                  int             n_comp,
                                  const Geometry& geom)
    :
    m_ngrow(n_grow),
    m_empty(false),
    m_initialized(false)
{
    initialize(ba,n_grow,n_comp,geom);
}

void
AuxBoundaryData::copy (const AuxBoundaryData& src,
                       int                    src_comp,
                       int                    dst_comp,
                       int                    num_comp)
{
    if (m_empty || src.m_empty) return;

    BL_ASSERT(m_initialized);
    BL_ASSERT(src_comp + num_comp <= src.m_fabs.nComp());
    BL_ASSERT(dst_comp + num_comp <= m_fabs.nComp());

    m_fabs.copy(src.m_fabs,src_comp,dst_comp,num_comp);
}

AuxBoundaryData::AuxBoundaryData (const AuxBoundaryData& rhs)
    :
    m_fabs(rhs.m_fabs.boxArray(),rhs.m_fabs.DistributionMap(),rhs.m_fabs.nComp(),0,
           MFInfo(), FArrayBoxFactory()),
    m_ngrow(rhs.m_ngrow)
{
    m_fabs.copy(rhs.m_fabs,0,0,rhs.m_fabs.nComp());

    m_initialized = true;
}

void
AuxBoundaryData::initialize (const BoxArray& ba,
			     int             n_grow,
			     int             n_comp,
                             const Geometry& geom)
{
    BL_ASSERT(!m_initialized);

    const bool verbose   = false;
    const int  NProcs    = ParallelDescriptor::NProcs();
    const Real strt_time = amrex::second();

    m_ngrow = n_grow;

    BoxList gcells = amrex::GetBndryCells(ba,n_grow);
    //
    // Remove any intersections with periodically shifted valid region.
    //
    if (geom.isAnyPeriodic())
    {
        Box dmn = geom.Domain();

        for (int d = 0; d < AMREX_SPACEDIM; d++) {
            if (!geom.isPeriodic(d)) {
                dmn.grow(d,n_grow);
            }
        }

        gcells.intersect(dmn);
    }

    gcells.simplify();

    if (gcells.size() < NProcs)
    {
        gcells.maxSize(AMREX_SPACEDIM == 3 ? 64 : 128);
    }

    BoxArray nba(gcells);
    DistributionMapping ndm{nba};

    gcells.clear();

    if (nba.size() > 0)
    {
        m_fabs.define(nba, ndm, n_comp, 0, MFInfo(), FArrayBoxFactory());
    }
    else
    {
        m_empty = true;
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = amrex::second() - strt_time;
	const int sz       = nba.size();

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);
	amrex::Print() << "AuxBoundaryData::initialize() size = " << sz << ", time = " << run_time << '\n';
#ifdef BL_LAZY
	});
#endif
    }

    m_initialized = true;
}

void
AuxBoundaryData::copyTo (MultiFab& mf,
                         int       src_comp,
                         int       dst_comp,
                         int       num_comp) const
{
    BL_ASSERT(m_initialized);

    if (!m_empty && mf.size() > 0)
    {
        mf.copy(m_fabs,src_comp,dst_comp,num_comp,0,mf.nGrow());
    }
}

void
AuxBoundaryData::copyFrom (const MultiFab& mf,
                           int       src_comp,
                           int       dst_comp,
                           int       num_comp,
			   int       src_ng)
{
    BL_ASSERT(m_initialized);

    if (!m_empty && mf.size() > 0)
    {
	m_fabs.copy(mf,src_comp,dst_comp,num_comp,src_ng,0);
    }
}

}
