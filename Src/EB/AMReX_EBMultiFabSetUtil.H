#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiFabSet.H>

namespace amrex
{
    void EB_set_covered (MultiFabSet& mfs,                                               Real   val) {
        for (int i = 0; i < mfs.nSet(); ++i) {
            EB_set_covered(mfs[i], val);
        }
    }
    void EB_set_covered (MultiFabSet& mfs, int icomp, int ncomp, int ngrow,              Real   val) {
        for (int i = 0; i < mfs.nSet(); ++i) {
            EB_set_covered(mfs[i], icomp, ncomp, ngrow, val);
        }
    }
    void EB_set_covered (MultiFabSet& mfs, int icomp, int ncomp,            const Vector<Real>& vals) {
        for (int i = 0; i < mfs.nSet(); ++i) {
            EB_set_covered(mfs[i], icomp, ncomp, vals);
        }
    }
    void EB_set_covered (MultiFabSet& mfs, int icomp, int ncomp, int ngrow, const Vector<Real>& vals) {
        for (int i = 0; i < mfs.nSet(); ++i) {
            EB_set_covered(mfs[i], icomp, ncomp, ngrow, vals);
        }
    }

    // void EB_set_covered_faces (const Array<MultiFabSet*,AMREX_SPACEDIM>& umac, Real val);
    // void EB_set_covered_faces (const Array<MultiFabSet*,AMREX_SPACEDIM>& umac, int scomp, int ncomp, const Vector<Real>& vals);

    void EB_average_down (const MultiFabSet& S_fine, MultiFabSet& S_crse, const MultiFabSet& vol_fine,
                          const MultiFabSet& vfrac_fine, int scomp, int ncomp, const IntVect& ratio) {
        for (int i = 0; i < S_fine.nSet(); ++i) {
            EB_average_down(S_fine[i], S_crse[i], vol_fine[i], vfrac_fine[i], scomp, ncomp, ratio);
        }
    }

    void EB_average_down (const MultiFabSet& S_fine, MultiFabSet& S_crse, int scomp, int ncomp,
                          int ratio) {
        for (int i = 0; i < S_fine.nSet(); ++i) {
            EB_average_down(S_fine[i], S_crse[i], scomp, ncomp, ratio);
        }
    }
    void EB_average_down (const MultiFabSet& S_fine, MultiFabSet& S_crse, int scomp, int ncomp,
                          const IntVect& ratio) {
        for (int i = 0; i < S_fine.nSet(); ++i) {
            EB_average_down(S_fine[i], S_crse[i], scomp, ncomp, ratio);
        }
    }

    // void EB_average_down_faces (const Array<const MultiFabSet*,AMREX_SPACEDIM>& fine,
    //                             const Array<      MultiFabSet*,AMREX_SPACEDIM>& crse,
    //                             int ratio, int ngcrse);
    // void EB_average_down_faces (const Array<const MultiFabSet*,AMREX_SPACEDIM>& fine,
    //                             const Array<      MultiFabSet*,AMREX_SPACEDIM>& crse,
    //                             const IntVect& ratio, int ngcrse);

    // //  This version takes periodicity into account.
    // void EB_average_down_faces (const Array<const MultiFabSet*,AMREX_SPACEDIM>& fine,
    //                             const Array<      MultiFabSet*,AMREX_SPACEDIM>& crse,
    //                             const IntVect& ratio, const Geometry& crse_geom);

    // void EB_average_down_boundaries (const MultiFabSet& fine, MultiFabSet& crse,
    //                                  int ratio, int ngcrse);
    // void EB_average_down_boundaries (const MultiFabSet& fine, MultiFabSet& crse,
    //                                  const IntVect& ratio, int ngcrse);

    // void EB_computeDivergence (MultiFabSet& divu, const Array<MultiFabSet const*,AMREX_SPACEDIM>& umac,
    //                            const Geometry& geom, bool already_on_centroids);

    // void EB_computeDivergence (MultiFabSet& divu, const Array<MultiFabSet const*,AMREX_SPACEDIM>& umac,
    //                            const Geometry& geom, bool already_on_centroids,
    //                            const MultiFabSet& vel_eb);

    // // Cell faces to cell centers
    // void EB_average_face_to_cellcenter (MultiFabSet& ccmf, int dcomp,
    //                                     const Array<MultiFabSet const*,AMREX_SPACEDIM>& fmf);

    // // Cell centers to cell centroids
    // void EB_interp_CC_to_Centroid (MultiFabSet& cent, const MultiFabSet& cc, int scomp, int dcomp, int ncomp, const Geometry& geom);

    // // Cell centers to face centroids
    // void EB_interp_CC_to_FaceCentroid (const MultiFabSet& cc,
    //                                    AMREX_D_DECL( MultiFabSet& fc_x,
    //                                                  MultiFabSet& fc_y,
    //                                                  MultiFabSet& fc_z),
    //                                    int scomp, int dcomp, int nc,
    //                                    const Geometry& geom,
    //                                    const amrex::Vector<amrex::BCRec>& a_bcs);

    // // Cell centroids to face centroids
    // void EB_interp_CellCentroid_to_FaceCentroid (const MultiFabSet& phi_centroid,
    //                                              const Array<MultiFabSet *,AMREX_SPACEDIM>& phi_faces,
    //                                              int scomp, int dcomp, int nc,
    //                                              const Geometry& geom,
    //                                              const amrex::Vector<amrex::BCRec>& a_bcs);

    // // Cell centroids to face centroids
    // void EB_interp_CellCentroid_to_FaceCentroid (const MultiFabSet& phi_centroid,
    //                                              const Vector<MultiFabSet*>& phi_faces,
    //                                              int scomp, int dcomp, int nc,
    //                                              const Geometry& geom,
    //                                              const amrex::Vector<amrex::BCRec>& a_bcs);

    // // Cell centroids to face centroids
    // void EB_interp_CellCentroid_to_FaceCentroid (const MultiFabSet& phi_centroid,
    //                                              AMREX_D_DECL( MultiFabSet& phi_xface,
    //                                                            MultiFabSet& phi_yface,
    //                                                            MultiFabSet& phi_zface),
    //                                              int scomp, int dcomp, int nc,
    //                                              const Geometry& geom,
    //                                              const amrex::Vector<amrex::BCRec>& a_bcs);

}

namespace amrex {

void
EB_set_covered (MultiFabSet& mf, Real val) {
    for (int s = 0; s < mf.nSet(); ++s) {
        amrex::EB_set_covered(mf[s], val);
    }
}

void
EB_set_covered (MultiFabSet& mf, int icomp, int ncomp, int ngrow,              Real   val);

void
EB_set_covered (MultiFabSet& mf, int icomp, int ncomp,            const Vector<Real>& vals);

void
EB_set_covered (MultiFabSet& mf, int icomp, int ncomp, int ngrow, const Vector<Real>& vals);

void
EB_set_covered_faces (const Array<MultiFabSet*,AMREX_SPACEDIM>& umac, Real val);

void
EB_set_covered_faces (const Array<MultiFabSet*,AMREX_SPACEDIM>& umac, int scomp, int ncomp, const Vector<Real>& vals)
{
    
}

}