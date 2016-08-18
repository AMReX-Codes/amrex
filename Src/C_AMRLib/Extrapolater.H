#ifndef _EXTRAPOLATER_H_
#define _EXTRAPOLATER_H_

#include <MultiFab.H>
#include <iMultiFab.H>

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

    void FirstOrderExtrap (MultiFab& mf, int scomp, int ncomp, const iMultiFab& mask);
}

#endif
