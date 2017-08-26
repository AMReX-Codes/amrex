#ifndef AMREX_BCUTIL_H_
#define AMREX_BCUTIL_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>

namespace amrex
{
    void FillDomainBoundary (MultiFab& phi, const Geometry& geom, const Array<BCRec>& bc);
}

#endif
