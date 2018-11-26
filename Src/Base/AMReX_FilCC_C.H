#ifndef AMREX_FILL_CC_H_
#define AMREX_FILL_CC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>

namespace amrex {

    AMREX_GPU_DEVICE
    void filcc_cell (const IntVect& iv, FArrayBox& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp);

}

#endif
