#ifndef AMREX_NONLOCAL_BC_H_
#define AMREX_NONLOCAL_BC_H_
#include <AMReX_Config.H>

#include <AMReX_FabArray.H>

namespace amrex { namespace NonLocalBC {

// Rotate90 fills the lo-x and lo-y boundary regions by rotating the data
// around (x=0,y=0) by 90 degrees in either direction.  It also fills the
// corner of lo-x and lo-y boundary region by rotating the data by 180
// degrees.

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
Rotate90 (FabArray<FAB>& mf, int scomp, int ncomp, IntVect const& nghost, Box const& domain);

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
Rotate90 (FabArray<FAB>& mf, Box const& domain);

// Rotate180 fills the lo-x boundary by rotating the data around
// (x=0,y=L_y/2) by 180 degrees.

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
Rotate180 (FabArray<FAB>& mf, int scomp, int ncomp, IntVect const& nghost, Box const& domain);

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
Rotate180 (FabArray<FAB>& mf, Box const& domain);

// Fill the polar boundaries of the spherical coordinates (theta, phi, r).
// The lo-x boundary is filled with f(-x,y) = f(x,mod(y+pi,2*pi)), and
// the hi-x boundary is filled with f(pi+x,y) = f(pi-x,mod(y+pi,2*pi)).

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
FillPolar (FabArray<FAB>& mf, int scomp, int ncomp, IntVect const& nghost, Box const& domain);

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
FillPolar (FabArray<FAB>& mf, Box const& domain);

}}

#include <AMReX_NonLocalBCImpl.H>

#endif
