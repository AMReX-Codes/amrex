
#include <winstd.H>

#include <iostream>

#include <AuxBoundaryData.H>

AuxBoundaryData::AuxBoundaryData ()
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
    m_fabs(rhs.m_fabs.boxArray(),rhs.m_fabs.nComp(),0),
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

    const bool verbose   = true;
    const int  NProcs    = ParallelDescriptor::NProcs();
    const Real strt_time = ParallelDescriptor::second();

    m_ngrow = n_grow;

    BoxList gcells = BoxLib::GetBndryCells(ba,n_grow);
    //
    // Remove any intersections with periodically shifted valid region.
    //
    if (geom.isAnyPeriodic())
    {
        Box dmn = geom.Domain();

        for (int d = 0; d < BL_SPACEDIM; d++)
            if (!geom.isPeriodic(d)) 
                dmn.grow(d,n_grow);

        for (BoxList::iterator it = gcells.begin(); it != gcells.end(); )
        {
            const Box isect = *it & dmn;

            if (isect.ok())
            {
                *it++ = isect;
            }
            else
            {
                gcells.remove(it++);
            }
        }
    }

    gcells.simplify();

    if (gcells.size() < NProcs)
    {
        gcells.maxSize(BL_SPACEDIM == 3 ? 64 : 128);
    }

    BoxArray nba(gcells);

    gcells.clear();

    if (nba.size() > 0)
    {
        std::vector<long> wgts;

        wgts.reserve(nba.size());

        for (int i = 0; i < nba.size(); i++)
            wgts.push_back(nba[i].numPts());

        DistributionMapping dm;
        //
        // This call doesn't invoke the MinimizeCommCosts() stuff.
        // There's very little to gain with this type of covering.
        // This also guarantees that this DM won't be put into the cache.
        //
        dm.KnapSackProcessorMap(wgts,NProcs);

        m_fabs.define(nba, n_comp, 0, dm, Fab_allocate);
    }
    else
    {
        m_empty = true;
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor()) 
            std::cout << "AuxBoundaryData::initialize() size = " << nba.size() << ", time = " << run_time << '\n';
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
        mf.copy(m_fabs,src_comp,dst_comp,num_comp);
    }
}

void
AuxBoundaryData::copyFrom (const MultiFab& mf,
                           int       src_comp,
                           int       dst_comp,
                           int       num_comp)
{
    BL_ASSERT(m_initialized);

    if (!m_empty && mf.size() > 0)
    {
        m_fabs.copy(mf,src_comp,dst_comp,num_comp);
    }
}
