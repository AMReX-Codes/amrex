
#include <AMReX_AuxBoundaryData.H>
#include <AMReX_Print.H>

#ifdef BL_LAZY
#include <AMReX_Lazy.H>
#endif

namespace amrex {
// \cond CODEGEN

AuxBoundaryData::AuxBoundaryData (const BoxArray& ba,
                                  int             n_grow,
                                  int             n_comp,
                                  const Geometry& geom)
    :
    m_ngrow(n_grow)
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

    m_fabs.ParallelCopy(src.m_fabs,src_comp,dst_comp,num_comp);
}

AuxBoundaryData::AuxBoundaryData (const AuxBoundaryData& rhs)
    :
    m_fabs(rhs.m_fabs.boxArray(),rhs.m_fabs.DistributionMap(),rhs.m_fabs.nComp(),0,
           MFInfo(), FArrayBoxFactory()),
    m_ngrow(rhs.m_ngrow),
    m_initialized(true)
{
    m_fabs.ParallelCopy(rhs.m_fabs,0,0,rhs.m_fabs.nComp());
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
    const auto strt_time = amrex::second(); // NOLINT

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

    if (!nba.empty())
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
        auto      run_time = amrex::second() - strt_time;
        const int sz       = static_cast<int>(nba.size());

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

    if (!m_empty && !mf.empty())
    {
        mf.ParallelCopy(m_fabs,src_comp,dst_comp,num_comp,0,mf.nGrow());
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

    if (!m_empty && !mf.empty())
    {
        m_fabs.ParallelCopy(mf,src_comp,dst_comp,num_comp,src_ng,0);
    }
}
// \endcond
}
