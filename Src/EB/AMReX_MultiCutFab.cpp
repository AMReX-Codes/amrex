
#include <AMReX_MultiCutFab.H>
#include <AMReX_MultiFab.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

std::size_t
CutFab::copyFromMem (const Box&  dstbox,
                     int         dstcomp,
                     int         numcomp,
                     const void* src)
{
    if (dptr != nullptr) {
        return FArrayBox::copyFromMem(dstbox, dstcomp, numcomp, src);
    } else {
        return sizeof(CutFab::value_type)*static_cast<std::size_t>(dstbox.numPts()*numcomp);
    }
}


CutFab&
CutFab::copy (const CutFab & src,
              const Box&     srcbox,
              int            srccomp,
              const Box&     destbox,
              int            destcomp,
              int            numcomp)
{
    if (dptr != nullptr) {
        FArrayBox::copy(src,srcbox,srccomp,destbox,destcomp,numcomp);
    }
    return *this;
}

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
            m_data.setFab(mfi, new CutFab(), false);
        }
    }
}

const CutFab&
MultiCutFab::operator[] (const MFIter& mfi) const
{
    AMREX_ASSERT(ok(mfi));
    return m_data[mfi];
}

CutFab&
MultiCutFab::operator[] (const MFIter& mfi)
{
    AMREX_ASSERT(ok(mfi));
    return m_data[mfi];
}

bool
MultiCutFab::ok (const MFIter& mfi) const
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
MultiCutFab::ParallelCopy (const MultiCutFab& src, int scomp, int dcomp, int ncomp, int sng, int dng)
{
    m_data.ParallelCopy(src.m_data, scomp, dcomp, ncomp, sng, dng);
}

void
MultiCutFab::ParallelCopy (const MultiCutFab& src, int scomp, int dcomp, int ncomp, int sng, int dng, const Periodicity& period)
{
    m_data.ParallelCopy(src.m_data, scomp, dcomp, ncomp, sng, dng, period);
}

}
