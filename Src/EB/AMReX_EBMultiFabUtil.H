#ifndef AMREX_EBMULTIFABUTIL_H_
#define AMREX_EBMULTIFABUTIL_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>

namespace amrex
{
    void EB_set_covered (MultiFab& mf,                                               Real   val);    
    void EB_set_covered (MultiFab& mf, int icomp, int ncomp, int ngrow,              Real   val);
    void EB_set_covered (MultiFab& mf, int icomp, int ncomp,            const Vector<Real>& vals);
    void EB_set_covered (MultiFab& mf, int icomp, int ncomp, int ngrow, const Vector<Real>& vals);

    void EB_set_covered_faces (const Array<MultiFab*,AMREX_SPACEDIM>& umac, Real val);

    void EB_average_down (const MultiFab& S_fine, MultiFab& S_crse, const MultiFab& vol_fine,
                          const MultiFab& vfrac_fine, int scomp, int ncomp, const IntVect& ratio);

    void EB_average_down (const MultiFab& S_fine, MultiFab& S_crse, int scomp, int ncomp,
                          int ratio);
    void EB_average_down (const MultiFab& S_fine, MultiFab& S_crse, int scomp, int ncomp,
                          const IntVect& ratio);

    void EB_average_down_faces (const Array<const MultiFab*,AMREX_SPACEDIM>& fine,
                                const Array<MultiFab*,AMREX_SPACEDIM>& crse,
                                int ratio, int ngcrse);
    void EB_average_down_faces (const Array<const MultiFab*,AMREX_SPACEDIM>& fine,
                                const Array<MultiFab*,AMREX_SPACEDIM>& crse,
                                const IntVect& ratio, int ngcrse);

    void EB_average_down_boundaries (const MultiFab& fine, MultiFab& crse,
                                     int ratio, int ngcrse);
    void EB_average_down_boundaries (const MultiFab& fine, MultiFab& crse,
                                     const IntVect& ratio, int ngcrse);

    void EB_computeDivergence (MultiFab& divu, const Array<MultiFab const*,AMREX_SPACEDIM>& umac,
                               const Geometry& geom);

    void EB_average_face_to_cellcenter (MultiFab& ccmf, int dcomp,
                                        const Array<MultiFab const*,AMREX_SPACEDIM>& fmf);
}

#endif
