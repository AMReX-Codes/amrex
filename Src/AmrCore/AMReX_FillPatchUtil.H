#ifndef AMREX_FillPatchUtil_H_
#define AMREX_FillPatchUtil_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MFInterpolater.H>
#include <AMReX_Array.H>
#include <AMReX_Utility.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBInterpolater.H>
#include <AMReX_EBMFInterpolater.H>
#endif

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <cmath>
#include <limits>
#include <map>

namespace amrex
{

    template <typename MFFAB>
    struct NullInterpHook
    {
        template <class F=MFFAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
        void operator() (MFFAB& /*fab*/, const Box& /*bx*/, int /*icomp*/, int /*ncomp*/) const {}

        template <class F=MFFAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
        void operator() (Array<MFFAB*, AMREX_SPACEDIM> /*fab*/, const Box& /*bx*/, int /*icomp*/, int /*ncomp*/) const {}

        template <class F=MFFAB, std::enable_if_t<IsFabArray<F>::value,int> = 0>
        void operator() (MFFAB& /*mf*/, int /*icomp*/, int /*ncomp*/) const {}
    };

    /**
     * \brief Test if AMR grids are properly nested.
     *
     * If grids are not properly nested, FillPatch functions may fail.
     *
     * \tparam Interp Interpolater type
     *
     * \param ratio refinement ratio
     * \param blocking_factor blocking factor on the fine level
     * \param ngrow number of ghost cells of fine MultiFab
     * \param boxType index type
     * \param mapper an interpolater object
     */
    template <typename Interp>
    bool ProperlyNested (const IntVect& ratio, const IntVect& blocking_factor, int ngrow,
                         const IndexType& boxType, Interp* mapper);

    /**
     * \brief FillPatch with data from the current level
     *
     * The destination MultiFab/FabArray is on the same AMR level as the
     * source MultiFab/FabArray. Usually this can only be used on AMR level
     * 0, because filling fine level MF usually requires coarse level
     * data. If needed, interpolation in time is performed.
     *
     * \tparam MF the MultiFab/FabArray type
     * \tparam BC functor for filling physical boundaries
     *
     * \param mf destination MF
     * \param nghost number of ghost cells of mf needed to be filled
     * \param time time associated with mf
     * \param smf source MFs
     * \param stime times associated smf
     * \param scomp starting component of the source MFs
     * \param dcomp starting component of the destination MF
     * \param ncomp number of components
     * \param geom Geometry for this level
     * \param physbcf functor for physical boundaries
     * \param bcfcomp starting component for physbcf
     */
    template <typename MF, typename BC>
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchSingleLevel (MF& mf, IntVect const& nghost, Real time,
                          const Vector<MF*>& smf, const Vector<Real>& stime,
                          int scomp, int dcomp, int ncomp,
                          const Geometry& geom,
                          BC& physbcf, int bcfcomp);

    /**
     * \brief FillPatch with data from the current level
     *
     * The destination MultiFab/FabArray is on the same AMR level as the
     * source MultiFab/FabArray. Usually this can only be used on AMR level
     * 0, because filling fine level MF usually requires coarse level
     * data. If needed, interpolation in time is performed. Ghost cells of
     * the destination MF are not filled.
     *
     * \tparam MF the MultiFab/FabArray type
     * \tparam BC functor for filling physical boundaries
     *
     * \param mf destination MF
     * \param time time associated with mf
     * \param smf source MFs
     * \param stime times associated smf
     * \param scomp starting component of the source MFs
     * \param dcomp starting component of the destination MF
     * \param ncomp number of components
     * \param geom Geometry for this level
     * \param physbcf functor for physical boundaries
     * \param bcfcomp starting component for physbcf
     */
    template <typename MF, typename BC>
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchSingleLevel (MF& mf, Real time,
                          const Vector<MF*>& smf, const Vector<Real>& stime,
                          int scomp, int dcomp, int ncomp,
                          const Geometry& geom,
                          BC& physbcf, int bcfcomp);

    /**
     * \brief FillPatch with data from the current level and the level below.
     *
     * First, we fill the destination MultiFab/FabArray with the current
     * level data as much as possible. This may include interpolation in
     * time. For the rest of the destination MF, we fill them with the
     * coarse level data using interpolation in space (and in time if
     * needed).
     *
     * \tparam MF the MultiFab/FabArray type
     * \tparam BC functor for filling physical boundaries
     * \tparam Interp spatial interpolater
     * \tparam PreInterpHook pre-interpolation hook
     * \tparam PostInterpHook post-interpolation hook
     *
     * \param mf destination MF on the fine level
     * \param nghost number of ghost cells of mf needed to be filled
     * \param time time associated with mf
     * \param cmf source MFs on the coarse level
     * \param ct times associated cmf
     * \param fmf source MFs on the fine level
     * \param ft times associated fmf
     * \param scomp starting component of the source MFs
     * \param dcomp starting component of the destination MF
     * \param ncomp number of components
     * \param cgeom Geometry for the coarse level
     * \param fgeom Geometry for the fine level
     * \param cbc functor for physical boundaries on the coarse level
     * \param cbccomp starting component for cbc
     * \param fbc functor for physical boundaries on the fine level
     * \param fbccomp starting component for fbc
     * \param ratio refinement ratio
     * \param mapper spatial interpolater
     * \param bcs boundary types for each component. We need this because some interpolaters need it.
     * \param bcscomp starting component for bcs
     * \param pre_interp pre-interpolation hook
     * \param post_interp post-interpolation hook
     */
    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels (MF& mf, IntVect const& nghost, Real time,
                        const Vector<MF*>& cmf, const Vector<Real>& ct,
                        const Vector<MF*>& fmf, const Vector<Real>& ft,
                        int scomp, int dcomp, int ncomp,
                        const Geometry& cgeom, const Geometry& fgeom,
                        BC& cbc, int cbccomp,
                        BC& fbc, int fbccomp,
                        const IntVect& ratio,
                        Interp* mapper,
                        const Vector<BCRec>& bcs, int bcscomp,
                        const PreInterpHook& pre_interp = {},
                        const PostInterpHook& post_interp = {});

    /**
     * \brief FillPatch with data from the current level and the level below.
     *
     * First, we fill the destination MultiFab/FabArray with the current
     * level data as much as possible. This may include interpolation in
     * time. For the rest of the destination MF, we fill them with the
     * coarse level data using interpolation in space (and in time if
     * needed). Ghost cells of the destination MF are not filled.
     *
     * \tparam MF the MultiFab/FabArray type
     * \tparam BC functor for filling physical boundaries
     * \tparam Interp spatial interpolater
     * \tparam PreInterpHook pre-interpolation hook
     * \tparam PostInterpHook post-interpolation hook
     *
     * \param mf destination MF on the fine level
     * \param time time associated with mf
     * \param cmf source MFs on the coarse level
     * \param ct times associated cmf
     * \param fmf source MFs on the fine level
     * \param ft times associated fmf
     * \param scomp starting component of the source MFs
     * \param dcomp starting component of the destination MF
     * \param ncomp number of components
     * \param cgeom Geometry for the coarse level
     * \param fgeom Geometry for the fine level
     * \param cbc functor for physical boundaries on the coarse level
     * \param cbccomp starting component for cbc
     * \param fbc functor for physical boundaries on the fine level
     * \param fbccomp starting component for fbc
     * \param ratio refinement ratio
     * \param mapper spatial interpolater
     * \param bcs boundary types for each component. We need this because some interpolaters need it.
     * \param bcscomp starting component for bcs
     * \param pre_interp pre-interpolation hook
     * \param post_interp post-interpolation hook
     */
    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels (MF& mf, Real time,
                        const Vector<MF*>& cmf, const Vector<Real>& ct,
                        const Vector<MF*>& fmf, const Vector<Real>& ft,
                        int scomp, int dcomp, int ncomp,
                        const Geometry& cgeom, const Geometry& fgeom,
                        BC& cbc, int cbccomp,
                        BC& fbc, int fbccomp,
                        const IntVect& ratio,
                        Interp* mapper,
                        const Vector<BCRec>& bcs, int bcscomp,
                        const PreInterpHook& pre_interp = {},
                        const PostInterpHook& post_interp = {});

    /**
     * \brief FillPatch for face variables with data from the current level
     * and the level below. Sometimes, we need to fillpatch all
     * AMREX_SPACEDIM face MultiFabs togother to satisfy certain constraint
     * such as divergence preserving.
     *
     * First, we fill the destination MultiFab/FabArray's with the current
     * level data as much as possible. This may include interpolation in
     * time. For the rest of the destination MFs, we fill them with the
     * coarse level data using interpolation in space (and in time if
     * needed).
     *
     * \tparam MF the MultiFab/FabArray type
     * \tparam BC functor for filling physical boundaries
     * \tparam Interp spatial interpolater
     * \tparam PreInterpHook pre-interpolation hook
     * \tparam PostInterpHook post-interpolation hook
     *
     * \param mf destination MFs on the fine level
     * \param nghost number of ghost cells of mf needed to be filled
     * \param time time associated with mf
     * \param cmf source MFs on the coarse level
     * \param ct times associated cmf
     * \param fmf source MFs on the fine level
     * \param ft times associated fmf
     * \param scomp starting component of the source MFs
     * \param dcomp starting component of the destination MFs
     * \param ncomp number of components
     * \param cgeom Geometry for the coarse level
     * \param fgeom Geometry for the fine level
     * \param cbc functor for physical boundaries on the coarse level
     * \param cbccomp starting component for cbc
     * \param fbc functor for physical boundaries on the fine level
     * \param fbccomp starting component for fbc
     * \param ratio refinement ratio
     * \param mapper spatial interpolater
     * \param bcs boundary types for each component. We need this because some interpolaters need it.
     * \param bcscomp starting component for bcs
     * \param pre_interp pre-interpolation hook
     * \param post_interp post-interpolation hook
     */
    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels (Array<MF*, AMREX_SPACEDIM> const& mf, IntVect const& nghost, Real time,
                        const Vector<Array<MF*, AMREX_SPACEDIM> >& cmf, const Vector<Real>& ct,
                        const Vector<Array<MF*, AMREX_SPACEDIM> >& fmf, const Vector<Real>& ft,
                        int scomp, int dcomp, int ncomp,
                        const Geometry& cgeom, const Geometry& fgeom,
                        Array<BC, AMREX_SPACEDIM>& cbc, const Array<int, AMREX_SPACEDIM>& cbccomp,
                        Array<BC, AMREX_SPACEDIM>& fbc, const Array<int, AMREX_SPACEDIM>& fbccomp,
                        const IntVect& ratio,
                        Interp* mapper,
                        const Array<Vector<BCRec>, AMREX_SPACEDIM>& bcs, const Array<int, AMREX_SPACEDIM>& bcscomp,
                        const PreInterpHook& pre_interp = {},
                        const PostInterpHook& post_interp = {});

    /**
     * \brief FillPatch for face variables with data from the current level
     * and the level below. Sometimes, we need to fillpatch all
     * AMREX_SPACEDIM face MultiFabs togother to satisfy certain constraint
     * such as divergence preserving.
     *
     * First, we fill the destination MultiFab/FabArray's with the current
     * level data as much as possible. This may include interpolation in
     * time. For the rest of the destination MFs, we fill them with the
     * coarse level data using interpolation in space (and in time if
     * needed).
     *
     * \tparam MF the MultiFab/FabArray type
     * \tparam BC functor for filling physical boundaries
     * \tparam Interp spatial interpolater
     * \tparam PreInterpHook pre-interpolation hook
     * \tparam PostInterpHook post-interpolation hook
     *
     * \param mf destination MFs on the fine level
     * \param nghost number of ghost cells of mf needed to be filled
     * \param time time associated with mf
     * \param cmf source MFs on the coarse level
     * \param ct times associated cmf
     * \param fmf source MFs on the fine level
     * \param ft times associated fmf
     * \param scomp starting component of the source MFs
     * \param dcomp starting component of the destination MFs
     * \param ncomp number of components
     * \param cgeom Geometry for the coarse level
     * \param fgeom Geometry for the fine level
     * \param cbc functor for physical boundaries on the coarse level
     * \param cbccomp starting component for cbc
     * \param fbc functor for physical boundaries on the fine level
     * \param fbccomp starting component for fbc
     * \param ratio refinement ratio
     * \param mapper spatial interpolater
     * \param bcs boundary types for each component. We need this because some interpolaters need it.
     * \param bcscomp starting component for bcs
     * \param pre_interp pre-interpolation hook
     * \param post_interp post-interpolation hook
     */
    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels (Array<MF*, AMREX_SPACEDIM> const& mf, IntVect const& nghost, Real time,
                        const Vector<Array<MF*, AMREX_SPACEDIM> >& cmf, const Vector<Real>& ct,
                        const Vector<Array<MF*, AMREX_SPACEDIM> >& fmf, const Vector<Real>& ft,
                        int scomp, int dcomp, int ncomp,
                        const Geometry& cgeom, const Geometry& fgeom,
                        Array<BC, AMREX_SPACEDIM>& cbc, int cbccomp,
                        Array<BC, AMREX_SPACEDIM>& fbc, int fbccomp,
                        const IntVect& ratio,
                        Interp* mapper,
                        const Array<Vector<BCRec>, AMREX_SPACEDIM>& bcs, int bcscomp,
                        const PreInterpHook& pre_interp = {},
                        const PostInterpHook& post_interp = {});

    /**
     * \brief FillPatch for face variables with data from the current level
     * and the level below. Sometimes, we need to fillpatch all
     * AMREX_SPACEDIM face MultiFabs togother to satisfy certain constraint
     * such as divergence preserving.
     *
     * First, we fill the destination MultiFab/FabArray's with the current
     * level data as much as possible. This may include interpolation in
     * time. For the rest of the destination MFs, we fill them with the
     * coarse level data using interpolation in space (and in time if
     * needed). Ghost cells of the destination MFs are not filled.
     *
     * \tparam MF the MultiFab/FabArray type
     * \tparam BC functor for filling physical boundaries
     * \tparam Interp spatial interpolater
     * \tparam PreInterpHook pre-interpolation hook
     * \tparam PostInterpHook post-interpolation hook
     *
     * \param mf destination MFs on the fine level
     * \param time time associated with mf
     * \param cmf source MFs on the coarse level
     * \param ct times associated cmf
     * \param fmf source MFs on the fine level
     * \param ft times associated fmf
     * \param scomp starting component of the source MFs
     * \param dcomp starting component of the destination MFs
     * \param ncomp number of components
     * \param cgeom Geometry for the coarse level
     * \param fgeom Geometry for the fine level
     * \param cbc functor for physical boundaries on the coarse level
     * \param cbccomp starting component for cbc
     * \param fbc functor for physical boundaries on the fine level
     * \param fbccomp starting component for fbc
     * \param ratio refinement ratio
     * \param mapper spatial interpolater
     * \param bcs boundary types for each component. We need this because some interpolaters need it.
     * \param bcscomp starting component for bcs
     * \param pre_interp pre-interpolation hook
     * \param post_interp post-interpolation hook
     */
    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels (Array<MF*, AMREX_SPACEDIM> const& mf, Real time,
                        const Vector<Array<MF*, AMREX_SPACEDIM> >& cmf, const Vector<Real>& ct,
                        const Vector<Array<MF*, AMREX_SPACEDIM> >& fmf, const Vector<Real>& ft,
                        int scomp, int dcomp, int ncomp,
                        const Geometry& cgeom, const Geometry& fgeom,
                        Array<BC, AMREX_SPACEDIM>& cbc, int cbccomp,
                        Array<BC, AMREX_SPACEDIM>& fbc, int fbccomp,
                        const IntVect& ratio,
                        Interp* mapper,
                        const Array<Vector<BCRec>, AMREX_SPACEDIM>& bcs, int bcscomp,
                        const PreInterpHook& pre_interp = {},
                        const PostInterpHook& post_interp = {});

#ifdef AMREX_USE_EB
    /**
     * \brief FillPatch with data from the current level and the level below.
     *
     * First, we fill the destination MultiFab/FabArray with the current
     * level data as much as possible. This may include interpolation in
     * time. For the rest of the destination MF, we fill them with the
     * coarse level data using interpolation in space (and in time if
     * needed).
     *
     * \tparam MF the MultiFab/FabArray type
     * \tparam BC functor for filling physical boundaries
     * \tparam Interp spatial interpolater
     * \tparam PreInterpHook pre-interpolation hook
     * \tparam PostInterpHook post-interpolation hook
     *
     * \param mf destination MF on the fine level
     * \param nghost number of ghost cells of mf needed to be filled
     * \param time time associated with mf
     * \param index_space EB IndexSpace
     * \param cmf source MFs on the coarse level
     * \param ct times associated cmf
     * \param fmf source MFs on the fine level
     * \param ft times associated fmf
     * \param scomp starting component of the source MFs
     * \param dcomp starting component of the destination MF
     * \param ncomp number of components
     * \param cgeom Geometry for the coarse level
     * \param fgeom Geometry for the fine level
     * \param cbc functor for physical boundaries on the coarse level
     * \param cbccomp starting component for cbc
     * \param fbc functor for physical boundaries on the fine level
     * \param fbccomp starting component for fbc
     * \param ratio refinement ratio
     * \param mapper spatial interpolater
     * \param bcs boundary types for each component. We need this because some interpolaters need it.
     * \param bcscomp starting component for bcs
     * \param pre_interp pre-interpolation hook
     * \param post_interp post-interpolation hook
     */
    template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels (MF& mf, IntVect const& nghost, Real time,
                        const EB2::IndexSpace& index_space,
                        const Vector<MF*>& cmf, const Vector<Real>& ct,
                        const Vector<MF*>& fmf, const Vector<Real>& ft,
                        int scomp, int dcomp, int ncomp,
                        const Geometry& cgeom, const Geometry& fgeom,
                        BC& cbc, int cbccomp,
                        BC& fbc, int fbccomp,
                        const IntVect& ratio,
                        Interp* mapper,
                        const Vector<BCRec>& bcs, int bcscomp,
                        const PreInterpHook& pre_interp,
                        const PostInterpHook& post_interp);

    /**
     * \brief FillPatch with data from the current level and the level below.
     *
     * First, we fill the destination MultiFab/FabArray with the current
     * level data as much as possible. This may include interpolation in
     * time. For the rest of the destination MF, we fill them with the
     * coarse level data using interpolation in space (and in time if
     * needed). Ghost cells of the destination MF are not filled.
     *
     * \tparam MF the MultiFab/FabArray type
     * \tparam BC functor for filling physical boundaries
     * \tparam Interp spatial interpolater
     * \tparam PreInterpHook pre-interpolation hook
     * \tparam PostInterpHook post-interpolation hook
     *
     * \param mf destination MF on the fine level
     * \param time time associated with mf
     * \param index_space EB IndexSpace
     * \param cmf source MFs on the coarse level
     * \param ct times associated cmf
     * \param fmf source MFs on the fine level
     * \param ft times associated fmf
     * \param scomp starting component of the source MFs
     * \param dcomp starting component of the destination MF
     * \param ncomp number of components
     * \param cgeom Geometry for the coarse level
     * \param fgeom Geometry for the fine level
     * \param cbc functor for physical boundaries on the coarse level
     * \param cbccomp starting component for cbc
     * \param fbc functor for physical boundaries on the fine level
     * \param fbccomp starting component for fbc
     * \param ratio refinement ratio
     * \param mapper spatial interpolater
     * \param bcs boundary types for each component. We need this because some interpolaters need it.
     * \param bcscomp starting component for bcs
     * \param pre_interp pre-interpolation hook
     * \param post_interp post-interpolation hook
     */
    template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels (MF& mf, Real time,
                        const EB2::IndexSpace& index_space,
                        const Vector<MF*>& cmf, const Vector<Real>& ct,
                        const Vector<MF*>& fmf, const Vector<Real>& ft,
                        int scomp, int dcomp, int ncomp,
                        const Geometry& cgeom, const Geometry& fgeom,
                        BC& cbc, int cbccomp,
                        BC& fbc, int fbccomp,
                        const IntVect& ratio,
                        Interp* mapper,
                        const Vector<BCRec>& bcs, int bcscomp,
                        const PreInterpHook& pre_interp,
                        const PostInterpHook& post_interp);
#endif

    /**
     * \brief Fill with interpolation of coarse level data
     *
     * No ghost cells of the destination MF are filled.
     *
     * \tparam MF the MultiFab/FabArray type
     * \tparam BC functor for filling physical boundaries
     * \tparam Interp spatial interpolater
     * \tparam PreInterpHook pre-interpolation hook
     * \tparam PostInterpHook post-interpolation hook
     *
     * \param mf destination MF on the fine level
     * \param time time associated with mf
     * \param cmf source MF on the coarse level
     * \param scomp starting component of the source MF
     * \param dcomp starting component of the destination MF
     * \param ncomp number of components
     * \param cgeom Geometry for the coarse level
     * \param fgeom Geometry for the fine level
     * \param cbc functor for physical boundaries on the coarse level
     * \param cbccomp starting component for cbc
     * \param fbc functor for physical boundaries on the fine level
     * \param fbccomp starting component for fbc
     * \param ratio refinement ratio
     * \param mapper spatial interpolater
     * \param bcs boundary types for each component. We need this because some interpolaters need it.
     * \param bcscomp starting component for bcs
     * \param pre_interp pre-interpolation hook
     * \param post_interp post-interpolation hook
     */
    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    InterpFromCoarseLevel (MF& mf, Real time,
                           const MF& cmf, int scomp, int dcomp, int ncomp,
                           const Geometry& cgeom, const Geometry& fgeom,
                           BC& cbc, int cbccomp,
                           BC& fbc, int fbccomp,
                           const IntVect& ratio,
                           Interp* mapper,
                           const Vector<BCRec>& bcs, int bcscomp,
                           const PreInterpHook& pre_interp = {},
                           const PostInterpHook& post_interp = {});

    /**
     * \brief Fill with interpolation of coarse level data
     *
     * \tparam MF the MultiFab/FabArray type
     * \tparam BC functor for filling physical boundaries
     * \tparam Interp spatial interpolater
     * \tparam PreInterpHook pre-interpolation hook
     * \tparam PostInterpHook post-interpolation hook
     *
     * \param mf destination MF on the fine level
     * \param nghost number of ghost cells of mf needed to be filled
     * \param time time associated with mf
     * \param cmf source MF on the coarse level
     * \param scomp starting component of the source MF
     * \param dcomp starting component of the destination MF
     * \param ncomp number of components
     * \param cgeom Geometry for the coarse level
     * \param fgeom Geometry for the fine level
     * \param cbc functor for physical boundaries on the coarse level
     * \param cbccomp starting component for cbc
     * \param fbc functor for physical boundaries on the fine level
     * \param fbccomp starting component for fbc
     * \param ratio refinement ratio
     * \param mapper spatial interpolater
     * \param bcs boundary types for each component. We need this because some interpolaters need it.
     * \param bcscomp starting component for bcs
     * \param pre_interp pre-interpolation hook
     * \param post_interp post-interpolation hook
     */
    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    InterpFromCoarseLevel (MF& mf, IntVect const& nghost, Real time,
                           const MF& cmf, int scomp, int dcomp, int ncomp,
                           const Geometry& cgeom, const Geometry& fgeom,
                           BC& cbc, int cbccomp,
                           BC& fbc, int fbccomp,
                           const IntVect& ratio,
                           Interp* mapper,
                           const Vector<BCRec>& bcs, int bcscomp,
                           const PreInterpHook& pre_interp = {},
                           const PostInterpHook& post_interp = {});

    /**
     * \brief Fill face variables with data from the coarse level.
     * Sometimes, we need to fillpatch all AMREX_SPACEDIM face MultiFabs
     * togother to satisfy certain constraint such as divergence preserving.
     *
     * \tparam MF the MultiFab/FabArray type
     * \tparam BC functor for filling physical boundaries
     * \tparam Interp spatial interpolater
     * \tparam PreInterpHook pre-interpolation hook
     * \tparam PostInterpHook post-interpolation hook
     *
     * \param mf destination MFs on the fine level
     * \param time time associated with mf
     * \param cmf source MFs on the coarse level
     * \param scomp starting component of the source MFs
     * \param dcomp starting component of the destination MFs
     * \param ncomp number of components
     * \param cgeom Geometry for the coarse level
     * \param fgeom Geometry for the fine level
     * \param cbc functor for physical boundaries on the coarse level
     * \param cbccomp starting component for cbc
     * \param fbc functor for physical boundaries on the fine level
     * \param fbccomp starting component for fbc
     * \param ratio refinement ratio
     * \param mapper spatial interpolater
     * \param bcs boundary types for each component. We need this because some interpolaters need it.
     * \param bcscomp starting component for bcs
     * \param pre_interp pre-interpolation hook
     * \param post_interp post-interpolation hook
     */
    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    InterpFromCoarseLevel (Array<MF*, AMREX_SPACEDIM> const& mf, Real time,
                           const Array<MF*, AMREX_SPACEDIM>& cmf, int scomp, int dcomp, int ncomp,
                           const Geometry& cgeom, const Geometry& fgeom,
                           Array<BC, AMREX_SPACEDIM>& cbc, int cbccomp,
                           Array<BC, AMREX_SPACEDIM>& fbc, int fbccomp,
                           const IntVect& ratio,
                           Interp* mapper,
                           const Array<Vector<BCRec>, AMREX_SPACEDIM>& bcs, int bcscomp,
                           const PreInterpHook& pre_interp = {},
                           const PostInterpHook& post_interp = {});

    /**
     * \brief Fill face variables with data from the coarse level.
     * Sometimes, we need to fillpatch all AMREX_SPACEDIM face MultiFabs
     * togother to satisfy certain constraint such as divergence preserving.
     *
     * \tparam MF the MultiFab/FabArray type
     * \tparam BC functor for filling physical boundaries
     * \tparam Interp spatial interpolater
     * \tparam PreInterpHook pre-interpolation hook
     * \tparam PostInterpHook post-interpolation hook
     *
     * \param mf destination MFs on the fine level
     * \param nghost number of ghost cells of mf needed to be filled
     * \param time time associated with mf
     * \param cmf source MFs on the coarse level
     * \param scomp starting component of the source MFs
     * \param dcomp starting component of the destination MFs
     * \param ncomp number of components
     * \param cgeom Geometry for the coarse level
     * \param fgeom Geometry for the fine level
     * \param cbc functor for physical boundaries on the coarse level
     * \param cbccomp starting component for cbc
     * \param fbc functor for physical boundaries on the fine level
     * \param fbccomp starting component for fbc
     * \param ratio refinement ratio
     * \param mapper spatial interpolater
     * \param bcs boundary types for each component. We need this because some interpolaters need it.
     * \param bcscomp starting component for bcs
     * \param pre_interp pre-interpolation hook
     * \param post_interp post-interpolation hook
     */
    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    InterpFromCoarseLevel (Array<MF*, AMREX_SPACEDIM> const& mf, IntVect const& nghost, Real time,
                           const Array<MF*, AMREX_SPACEDIM>& cmf, int scomp, int dcomp, int ncomp,
                           const Geometry& cgeom, const Geometry& fgeom,
                           Array<BC, AMREX_SPACEDIM>& cbc, int cbccomp,
                           Array<BC, AMREX_SPACEDIM>& fbc, int fbccomp,
                           const IntVect& ratio,
                           Interp* mapper,
                           const Array<Vector<BCRec>, AMREX_SPACEDIM>& bcs, int bcscomp,
                           const PreInterpHook& pre_interp = {},
                           const PostInterpHook& post_interp = {});

#ifndef BL_NO_FORT
    enum InterpEM_t { InterpE, InterpB};

    void InterpCrseFineBndryEMfield (InterpEM_t interp_type,
                                     const Array<MultiFab,AMREX_SPACEDIM>& crse,
                                     Array<MultiFab,AMREX_SPACEDIM>& fine,
                                     const Geometry& cgeom, const Geometry& fgeom,
                                     int ref_ratio);

    void InterpCrseFineBndryEMfield (InterpEM_t interp_type,
                                     const Array<MultiFab const*,AMREX_SPACEDIM>& crse,
                                     const Array<MultiFab*,AMREX_SPACEDIM>& fine,
                                     const Geometry& cgeom, const Geometry& fgeom,
                                     int ref_ratio);
#endif

    /**
     * \brief FillPatch with data from AMR levels.
     *
     * First, we try to fill the destination MultiFab/FabArray with this
     * level's data if it's available. For the unfilled region, we try to
     * fill with the coarse level below if it's available. Even coarser
     * levels will be used if necessary till all regions are filled. This
     * function is more expensive than FillPatchTwoLevels. So if one knows
     * FillPatchTwoLevels can do the job because grids are properly nested,
     * this function should be avoided.
     *
     * \tparam MF the MultiFab/FabArray type
     * \tparam BC functor for filling physical boundaries
     * \tparam Interp spatial interpolater
     *
     * \param mf destination MF
     * \param level AMR level associated with mf
     * \param nghost number of ghost cells of mf needed to be filled
     * \param time time associated with mf
     * \param smf source MFs. The outer Vector is for AMR levels, whereas the inner
     *            Vector is for data at various times. It is not an error if the
     *            level for the destination MF is finer than data in smf (i.e.,
     *            `level >= smf.size()`).
     * \param st times associated smf
     * \param scomp starting component of the source MFs
     * \param dcomp starting component of the destination MF
     * \param ncomp number of components
     * \param geom Geometry objects for AMR levels. The size must be big enough
     *             such that `level < geom.size()`.
     * \param bc functors for physical boundaries on AMR levels. The size must be
     *           big enough such that `level < bc.size()`.
     * \param bccomp starting component for bc
     * \param ratio refinement ratio for AMR levels. The size must be big enough
     *              such that `level < bc.size()-1`.
     * \param mapper spatial interpolater
     * \param bcr boundary types for each component. We need this because some interpolaters need it.
     * \param bcrcomp starting component for bcr
     */
    template <typename MF, typename BC, typename Interp>
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchNLevels (MF& mf, int level, const IntVect& nghost, Real time,
                      const Vector<Vector<MF*>>& smf, const Vector<Vector<Real>>& st,
                      int scomp, int dcomp, int ncomp,
                      const Vector<Geometry>& geom,
                      Vector<BC>& bc, int bccomp,
                      const Vector<IntVect>& ratio,
                      Interp* mapper,
                      const Vector<BCRec>& bcr, int bcrcomp);

}

#include <AMReX_FillPatchUtil_I.H>

#endif
