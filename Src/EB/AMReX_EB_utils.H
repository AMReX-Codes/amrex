#ifndef AMREX_EB_UTILS_H_
#define AMREX_EB_UTILS_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_GeometryShop.H>

namespace amrex {

#if (AMREX_SPACEDIM > 1)

    void single_level_redistribute (MultiFab& div_tmp_in, MultiFab& div_out,
                                    int div_comp, int ncomp, const Geometry& geom);

    void single_level_weighted_redistribute (MultiFab& div_tmp_in, MultiFab& div_out, const MultiFab& weights,
                                             int div_comp, int ncomp, const Geometry& geom);

#if (__cplusplus >= 201402L)
    [[deprecated]]
#endif
    inline void single_level_redistribute (int lev, MultiFab& div_tmp_in, MultiFab& div_out,
                                           int div_comp, int ncomp, const Vector<Geometry>& geom)
    { single_level_redistribute(div_tmp_in, div_out, div_comp, ncomp, geom[lev]); }

#if (__cplusplus >= 201402L)
    [[deprecated]]
#endif
    inline void single_level_weighted_redistribute (int lev, MultiFab& div_tmp_in, MultiFab& div_out, const MultiFab& weights,
                                                    int div_comp, int ncomp, const Vector<Geometry>& geom)
    { single_level_weighted_redistribute(div_tmp_in, div_out, weights, div_comp, ncomp, geom[lev]); }
#endif

    template <typename G>
    void FillImpFunc (MultiFab& mf, G const& gshop, Geometry const& geom)
    {
        AMREX_ALWAYS_ASSERT(mf.is_nodal());

        Box bounding_box = geom.Domain();
        bounding_box.surroundingNodes();
        bool extend_domain_face = EB2::ExtendDomainFace();
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (!extend_domain_face || geom.isPeriodic(idim)) {
                bounding_box.grow(mf.nGrowVect()[idim]);
            }
        }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for(MFIter mfi(mf); mfi.isValid(); ++ mfi) {
            gshop.fillFab(mf[mfi], geom, RunOn::Gpu, bounding_box);
        }
    }

#if (AMREX_SPACEDIM == 3)
    // mf is at the same level as ls_lev.
    // ls_lev is refratio finer than eb_factory.
    void FillSignedDistance (MultiFab& mf, EB2::Level const& ls_lev,
                             EBFArrayBoxFactory const& eb_fac, int refratio,
                             bool fluid_has_positive_sign=true);

    // mf needs to be built with EBFactory already.
    void FillSignedDistance (MultiFab& mf, bool fluid_has_positive_sign=true);
#endif
}

#endif
