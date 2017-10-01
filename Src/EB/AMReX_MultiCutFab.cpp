
#include <AMReX_MultiCutFab.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

MultiCutFab::MultiCutFab ()
{}

MultiCutFab::MultiCutFab (const BoxArray& ba, const DistributionMapping& dm,
                          int ncomp, int ngrow, const FabArray<EBCellFlagFab>& cellflags)
    : m_data(ba,dm,ncomp,ngrow),
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
    m_data.define(ba,dm,ncomp,ngrow);
    m_cellflags = &cellflags;
    remove();
}

void
MultiCutFab::remove ()
{
    for (MFIter mfi(m_data); mfi.isValid(); ++mfi)
    {
        if ((*m_cellflags)[mfi].getType() != FabType::singlevalued)
        {
            FArrayBox* p = &(m_data[mfi]);
            delete p;
            m_data.setFab(mfi, nullptr, false);
        }
    }
}

const FArrayBox&
MultiCutFab::operator[] (const MFIter& mfi) const
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

}
