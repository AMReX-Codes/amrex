#ifndef AMREX_EB_MULTIFABSET_UTIL_H_
#define AMREX_EB_MULTIFABSET_UTIL_H_

#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_FabArraySet.H>

namespace amrex
{
    void EB_set_covered (MultiFabSet& mfs,                                               Real   val);
    void EB_set_covered (MultiFabSet& mfs, int icomp, int ncomp, int ngrow,              Real   val);
    void EB_set_covered (MultiFabSet& mfs, int icomp, int ncomp,            const Vector<Real>& vals);
    void EB_set_covered (MultiFabSet& mfs, int icomp, int ncomp, int ngrow, const Vector<Real>& vals);

    void EB_average_down (const MultiFabSet& fine, MultiFabSet& crse,
                          const IntVect& ratio, const Vector<const IntVect>& ngrow);

    void EB_average_down_same_ixtype (const MultiFab& fine, MultiFab& crse, const IntVect& ratio,
                                      const IntVect& ngrow);

}


#endif /* AMREX_EB_MULTIFABSET_UTIL_H_ */