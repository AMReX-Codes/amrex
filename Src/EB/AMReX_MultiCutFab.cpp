
#include <AMReX_MultiCutFab.H>

namespace amrex {

MultiCutFab::MultiCutFab ()
{}

MultiCutFab::MultiCutFab (const BoxArray& ba, const DistributionMapping& dm,
                          int ncomp, int ngrow, const FabArray<EBCellFlagFab>& cellflags)
    : m_data(ba,dm,ncomp,ngrow),
      m_cellflags(&cellflags)
{}

MultiCutFab::~MultiCutFab ()
{}

void
MultiCutFab::define (const BoxArray& ba, const DistributionMapping& dm,
                     int ncomp, int ngrow, const FabArray<EBCellFlagFab>& cellflags)
{
    m_data.define(ba,dm,ncomp,ngrow);
    m_cellflags = &cellflags;
}

const FArrayBox&
MultiCutFab::operator[] (const MFIter& mfi) const
{
    AMREX_ASSERT((*m_cellflags)[mfi].getType() == FabType::singlevalued);
    return m_data[mfi];
}

}
