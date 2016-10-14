#ifndef _MultiFabUtil_H_
#define _MultiFabUtil_H_

#include <MultiFab.H>
#include <Geometry.H>

namespace BoxLib
{
    // Average face-based MultiFab onto cell-centered MultiFab.
    void average_edge_to_cellcenter (MultiFab& cc, int dcomp, const std::vector<MultiFab*>& edge);

    // Average face-based MultiFab onto cell-centered MultiFab.
    void average_face_to_cellcenter (MultiFab& cc, int dcomp, const std::vector<MultiFab*>& fc);
    void average_face_to_cellcenter (MultiFab& cc, const PArray<MultiFab>& fc, const Geometry& geom);

    // Average cell-centered MultiFab onto face-based MultiFab.
    void average_cellcenter_to_face (PArray<MultiFab>& fc, const MultiFab& cc, const Geometry& geom);

    // Average fine face-based MultiFab onto crse fine-centered MultiFab.
    // This routine assumes that the crse BoxArray is a coarsened version of the fine BoxArray.
    void average_down_faces (PArray<MultiFab>& fine, PArray<MultiFab>& crse, IntVect& ratio);

    // Average fine cell-based MultiFab onto crse cell-centered MultiFab using volume-weighting.
    // This routine DOES NOT assume that the crse BoxArray is a coarsened version of the fine BoxArray.
    void average_down(MultiFab& S_fine, MultiFab& S_crse, const Geometry& fine_geom, const Geometry& crse_geom,
                      int scomp, int ncomp, const IntVect& ratio);
    void average_down(MultiFab& S_fine, MultiFab& S_crse, const Geometry& fine_geom, const Geometry& crse_geom,
                      int scomp, int ncomp, const int ratio);

    // Average fine cell-based MultiFab onto crse cell-centered MultiFab without volume weighting.
    // This routine DOES NOT assume that the crse BoxArray is a coarsened version of the fine BoxArray.
    void average_down(MultiFab& S_fine, MultiFab& S_crse, int scomp, int ncomp, const IntVect& ratio);
    void average_down(MultiFab& S_fine, MultiFab& S_crse, int scomp, int ncomp,       int      ratio);

    // Fill MultiFab boundaries except for physical boundaries.
    void fill_boundary(MultiFab& mf, const Geometry& geom, bool cross = false);
    void fill_boundary(MultiFab& mf, int scomp, int ncomp, const Geometry& geom, bool cross = false);
}

#endif
