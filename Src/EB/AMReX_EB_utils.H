#ifndef AMREX_EB_UTILS_H_
#define AMREX_EB_UTILS_H_
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_GeometryShop.H>

namespace amrex {
    /**
     * \brief Fill MultiFab with implicit function.
     *
     * This function fills the nodal MultiFab with the implicit function in
     * GeometryShop.  Note that an implicit function is not necessarily a
     * signed distance function.
     *
     * \tparam G is the GeometryShop type
     *
     * \param mf is a nodal MultiFab.
     * \param gshop is a GeometryShop object.
     * \param geom is a Geometry object.
     */
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

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for(MFIter mfi(mf); mfi.isValid(); ++ mfi) {
            gshop.fillFab(mf[mfi], geom, RunOn::Gpu, bounding_box);
        }
    }

    /**
     * \brief Fill MultiFab with signed distance.
     *
     * This function fills the nodal MultiFab with signed distance.  Note
     * that the distance is valid only if it's within a few cells to the EB.
     *
     * \param mf is a nodal MultiFab.
     * \param ls_lev is an EB2::Level obejct with an implicit function. This
     *               is at the same level as mf.
     * \param eb_fac is an EBFArrayBoxFactory object containing EB informaiton.
     * \param refratio is the refinement ratio of mf to eb_fac.
     * \param fluid_has_positive_sign determines the sign of the fluid.
     */
    void FillSignedDistance (MultiFab& mf, EB2::Level const& ls_lev,
                             EBFArrayBoxFactory const& eb_fac, int refratio,
                             bool fluid_has_positive_sign=true);

    /**
     * \brief Fill MultiFab with signed distance.
     *
     * This function fills the nodal MultiFab with signed distance.  Note
     * that the distance is valid only if it's within a few cells to the EB.
     * The MultiFab must have been built with an EBFArrayBoxFactory.
     *
     * \param mf is a nodal MultiFab built with EBFArrayBoxFactory.
     * \param fluid_has_positive_sign determines the sign of the fluid.
     */
    void FillSignedDistance (MultiFab& mf, bool fluid_has_positive_sign=true);
}

#endif
