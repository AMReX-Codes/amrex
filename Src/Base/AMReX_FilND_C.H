// IWYU pragma: private, include "AMReX_PhysBCFunct.H"
#ifndef AMREX_FILL_ND_H_
#define AMREX_FILL_ND_H_
#include <AMReX_Config.H>

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_BCRec.H>

namespace amrex {
    void fab_filnd (Box const& bx, Array4<Real> const& q, int ncomp,
                    Box const& domain, Real const* dx, Real const* xlo,
                    BCRec const* bc);
}

#endif
