#include <AMReX_PlotfileData.H>
#include <AMReX_DataServices.H>

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
PlotfileData::probDomain (int level) const
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

MultiFab
PlotfileData::get (int level, const std::string& varname)
{
    MultiFab mf(m_ad->boxArray(level), m_ad->DistributionMap(level), 1, 0);
    m_ad->FillVar(mf, level, varname, 0);
    return mf;
}

MultiFab
PlotfileData::get (int level)
{
    MultiFab mf(boxArray(level), DistributionMap(level), nComp(), 0);
    const auto& varnames = varNames();
    for (int n = 0; n < nComp(); ++n) {
        m_ad->FillVar(mf, level, varnames[n], n);
    }
    return mf;
}

void
PlotfileData::fill (MultiFab& mf, int destcomp,
                    int level, const std::string& varname)
{
    m_ad->FillVar(mf, level, varname, destcomp);
}

void
PlotfileData::flush ()
{
    m_ad->FlushGrids();
}

void
PlotfileData::flush (const std::string& varname)
{
    m_ad->FlushGrids(m_ad->StateNumber(varname));
}

}
