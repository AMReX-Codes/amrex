#ifndef _EXTRAPOLATER_H_
#define _EXTRAPOLATER_H_

#include <MultiFab.H>
#include <Geometry.H>

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

    // It is expected that FillBoundary (w/ periodicity) has been called on mf.
    void FirstOrderExtrap (MultiFab& mf, const Geometry& geom, int scomp, int ncomp);
}

#endif
