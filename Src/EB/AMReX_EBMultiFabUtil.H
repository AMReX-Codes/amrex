#ifndef AMREX_EBMULTIFABUTIL_H_
#define AMREX_EBMULTIFABUTIL_H_

#include <AMReX_MultiFab.H>
#include <array>

namespace amrex
{
    void EB_set_covered (MultiFab& mf);    
    void EB_set_covered (MultiFab& mf, int icomp, int ncomp);
    void EB_set_covered (MultiFab& mf, int icomp, int ncomp, Real val);
    void EB_set_single_valued_cells (MultiFab& mf, int icomp, int ncomp, Real val);    

    void EB_average_down (const MultiFab& S_fine, MultiFab& S_crse, const MultiFab& vol_fine,
                          const MultiFab& vfrac_fine, int scomp, int ncomp, const IntVect& ratio);

}

#endif
