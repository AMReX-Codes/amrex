#ifndef AMREX_BCUTIL_H_
#define AMREX_BCUTIL_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>

namespace amrex
{
    // This is for filling cell-centered data outside physical domain
    // (excluding periodic boundaries).  It only fills
    // BCType::foextrap, BCType::hoextrap, BCType::reflect_even, and
    // BCType::reflect::odd.  It does not fill BCType::ext_dir (i.e.,
    // external Dirichlet).  If you have BCType::ext_dir, you can
    // copy, paste and rename this function, and then modify it.
    void FillDomainBoundary (MultiFab& phi, const Geometry& geom, const Vector<BCRec>& bc);
}

#endif
