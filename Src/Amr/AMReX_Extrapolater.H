#ifndef AMREX_EXTRAPOLATER_H_
#define AMREX_EXTRAPOLATER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

namespace amrex {

namespace Extrapolater
{
    // finebnd: boundary cells covered by fine cells (including periodically shifted fine cells)
    // crsebnd: boundary cells not covered by fine cells
    // physbnd: boundary cells outside the domain (excluding periodic boundaries)
    // interior: interior cells
    const int finebnd = 1;
    const int crsebnd = 0;
    const int physbnd = 0;
    const int interior = 1;

    //! It is expected that FillBoundary (w/ periodicity) has been called on mf.
    void FirstOrderExtrap (MultiFab& mf, const Geometry& geom, int scomp, int ncomp);
}

}

#endif
