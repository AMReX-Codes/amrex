#ifndef AMREX_FILL_ND_H_
#define AMREX_FILL_ND_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>

namespace amrex {
    void fab_filnd (Box const& bx, Array4<Real> const& q, int ncomp,
                    Box const& domain, Real const* dx, Real const* xlo,
                    BCRec const* bc);
}

#endif
