#ifndef _MultiFabUtil_H_
#define _MultiFabUtil_H_

#include <MultiFab.H>
#include <Geometry.H>

namespace BoxLib
{
    // Average face-based MultiFab onto cell-centered MultiFab.
    void average_face_to_cellcenter (MultiFab& cc, const PArray<MultiFab>& fc, const Geometry& geom);

    // Fill MultiFab boundaries except for physical boundaries.
    void fill_boundary(MultiFab& mf, const Geometry& geom, bool cross = false);
    void fill_boundary(MultiFab& mf, int scomp, int ncomp, const Geometry& geom, bool cross = false);
}

#endif
