#include <AMReX_PlotfileData.H>

namespace amrex {

namespace {
    bool initialized = false;
}

PlotfileData::PlotfileData (std::string const& plotfile_name)
{
    if (!initialized) {
        initialized = true;
        DataServices::SetBatchMode();
        AmrData::SetDimensionAgnostic(true);
        AmrData::SetVerbose(true);
    }

    m_ds.reset(new DataServices(plotfile_name, Amrvis::NEWPLT));
    AMREX_ALWAYS_ASSERT(m_ds->AmrDataOk());
    m_ad = &(m_ds->AmrDataRef());
}

PlotfileData::~PlotfileData ()
{}

int
PlotfileData::spaceDim () const
{
    return m_ad->SpaceDim();
}

Real
PlotfileData::time () const
{
    return m_ad->Time();
}

int
PlotfileData::finestLevel () const
{
    return m_ad->FinestLevel();
}

IntVect
PlotfileData::refRatio (int level) const
{
    return IntVect(m_ad->RefRatio()[level]);
}

const BoxArray&
PlotfileData::boxArray (int level) const
{
    return m_ad->boxArray(level);
}

const DistributionMapping&
PlotfileData::DistributionMap (int level) const
{
    return m_ad->DistributionMap(level);
}

int
PlotfileData::coordSys () const
{
    return m_ad->CoordSys();
}

Box
PlotfileData::ProbDomain (int level) const
{
    return m_ad->ProbDomain()[level];
}

Array<Real,AMREX_SPACEDIM>
PlotfileData::probSize () const
{
    const auto& d = m_ad->ProbSize();
    return {AMREX_D_DECL(d[0],d[1],d[2])};
}

Array<Real,AMREX_SPACEDIM>
PlotfileData::probLo () const
{
    const auto& d = m_ad->ProbLo();
    return {AMREX_D_DECL(d[0],d[1],d[2])};
}

Array<Real,AMREX_SPACEDIM>
PlotfileData::probHi () const
{
    const auto& d = m_ad->ProbHi();
    return {AMREX_D_DECL(d[0],d[1],d[2])};
}

Array<Real,AMREX_SPACEDIM>
PlotfileData::cellSize (int level) const
{
    const auto& d = m_ad->CellSize(level);
    return {AMREX_D_DECL(d[0],d[1],d[2])};
}

const Vector<std::string>&
PlotfileData::varNames () const
{
    return m_ad->PlotVarNames();
}

int
PlotfileData::nComp () const
{
    return m_ad->NComp();
}

IntVect
PlotfileData::nGrowVect () const
{
    return m_ad->NGrowVect();
}


}
