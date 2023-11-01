#include <AMReX_EBMultiFabSetUtil.H>

namespace amrex {

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

}