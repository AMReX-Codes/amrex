
#include <AMReX_MultiCutFab.H>

namespace amrex {

MultiCutFab::MultiCutFab ()
{}

MultiCutFab::MultiCutFab (const BoxArray& ba, const DistributionMapping& dm, int ncomp, int ngrow)
    : m_data(ba,dm,ncomp,ngrow)
{}

MultiCutFab::~MultiCutFab ()
{}

void
MultiCutFab::define (const BoxArray& ba, const DistributionMapping& dm, int ncomp, int ngrow)
{
    m_data.define(ba,dm,ncomp,ngrow);
}

const FArrayBox&
MultiCutFab::operator[] (const MFIter& mfi) const
{
    //
    return m_data[mfi];
}

}
