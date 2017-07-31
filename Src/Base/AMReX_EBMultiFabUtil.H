#ifndef AMREX_EBMULTIFABUTIL_H_
#define AMREX_EBMULTIFABUTIL_H_

#include <AMReX_MultiFab.H>

namespace amrex
{
    void EB_set_covered (MultiFab& mf);    
    void EB_set_covered (MultiFab& mf, int icomp, int ncomp);
}

#endif
