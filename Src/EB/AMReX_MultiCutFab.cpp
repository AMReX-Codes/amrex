
#include <AMReX_MultiCutFab.H>
#include <AMReX_MultiFab.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

MultiCutFab::MultiCutFab ()
{}

MultiCutFab::MultiCutFab (const BoxArray& ba, const DistributionMapping& dm,
                          int ncomp, int ngrow, const FabArray<EBCellFlagFab>& cellflags)
    : m_data(ba,dm,ncomp,ngrow,MFInfo(),DefaultFabFactory<CutFab>()),
      m_cellflags(&cellflags)
{
    remove();
}

MultiCutFab::~MultiCutFab ()
{}

void
MultiCutFab::define (const BoxArray& ba, const DistributionMapping& dm,
                     int ncomp, int ngrow, const FabArray<EBCellFlagFab>& cellflags)
{
    m_data.define(ba,dm,ncomp,ngrow,MFInfo(),DefaultFabFactory<CutFab>()),
    m_cellflags = &cellflags;
    remove();
}

void
MultiCutFab::remove ()
{
    for (MFIter mfi(m_data); mfi.isValid(); ++mfi)
    {
        if (!ok(mfi))
        {
            CutFab* p = &(m_data[mfi]);
            delete p;
            m_data.setFab(mfi, ::new CutFab(), false);
        }
    }
}

const CutFab&
MultiCutFab::operator[] (const MFIter& mfi) const noexcept
{
    AMREX_ASSERT(ok(mfi));
    return m_data[mfi];
}

CutFab&
MultiCutFab::operator[] (const MFIter& mfi) noexcept
{
    AMREX_ASSERT(ok(mfi));
    return m_data[mfi];
}

CutFab const*
MultiCutFab::fabPtr (const MFIter& mfi) const noexcept
{
    AMREX_ASSERT(ok(mfi));
    return m_data.fabPtr(mfi);
}

CutFab*
MultiCutFab::fabPtr (const MFIter& mfi) noexcept
{
    AMREX_ASSERT(ok(mfi));
    return m_data.fabPtr(mfi);
}

Array4<Real const>
MultiCutFab::array (const MFIter& mfi) const noexcept
{
    AMREX_ASSERT(ok(mfi));
    return m_data.array(mfi);
}

Array4<Real>
MultiCutFab::array (const MFIter& mfi) noexcept
{
    AMREX_ASSERT(ok(mfi));
    return m_data.array(mfi);
}

bool
MultiCutFab::ok (const MFIter& mfi) const noexcept
{
    return (*m_cellflags)[mfi].getType() == FabType::singlevalued;
}

void
MultiCutFab::setVal (Real val)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(m_data); mfi.isValid(); ++mfi)
    {
        if (ok(mfi)) {
            m_data[mfi].setVal(val);
        }
    }
}

void
MultiCutFab::ParallelCopy (const MultiCutFab& src, int scomp, int dcomp, int ncomp, int sng, int dng, const Periodicity& period)
{
    m_data.ParallelCopy(src.m_data, scomp, dcomp, ncomp, sng, dng, period);
}

MultiFab
MultiCutFab::ToMultiFab (Real regular_value, Real covered_value) const
{
    MultiFab mf(boxArray(), DistributionMap(), nComp(), nGrow());
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        auto t = (*m_cellflags)[mfi].getType();
        if (t == FabType::singlevalued) {
            mf[mfi].copy(m_data[mfi]);
        } else if (t == FabType::regular) {
            mf[mfi].setVal(regular_value);
        } else {
            mf[mfi].setVal(covered_value);
        }
    }
    return mf;
}

}
