#ifndef AMREX_FILL_CC_H_
#define AMREX_FILL_CC_H_

#if (AMREX_SPACEDIM == 1)
#include <AMReX_FilCC_1D_C.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_FilCC_2D_C.H>
#else
#include <AMReX_FilCC_3D_C.H>
#endif

namespace amrex {
    void fab_filcc (Box const& bx, Array4<Real> const& q, int ncomp,
                    Box const& domain, Real const* dx, Real const* xlo,
                    BCRec const* bc);
}

#endif
