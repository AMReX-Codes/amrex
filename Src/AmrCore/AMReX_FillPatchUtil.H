#ifndef AMREX_FillPatchUtil_H_
#define AMREX_FillPatchUtil_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Interpolater.H>
#include <AMReX_Array.H>
#include <AMReX_Utility.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>
#endif

#include <cmath>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex
{

    template <typename FAB>
    struct NullInterpHook
    {
        void operator() (FAB& /*fab*/, const Box& /*bx*/, int /*icomp*/, int /*ncomp*/) const {}
    };

    template <typename Interp>
    bool ProperlyNested (const IntVect& ratio, const IntVect& blocking_factor, int ngrow,
			 const IndexType& boxType, Interp* mapper);

    template <typename MF, typename BC>
    EnableIf_t<IsFabArray<MF>::value>
    FillPatchSingleLevel (MF& mf, IntVect const& nghost, Real time,
                          const Vector<MF*>& smf, const Vector<Real>& stime,
                          int scomp, int dcomp, int ncomp,
                          const Geometry& geom,
                          BC& physbcf, int bcfcomp);

    template <typename MF, typename BC>
    EnableIf_t<IsFabArray<MF>::value>
    FillPatchSingleLevel (MF& mf, Real time,
                          const Vector<MF*>& smf, const Vector<Real>& stime,
                          int scomp, int dcomp, int ncomp,
                          const Geometry& geom,
                          BC& physbcf, int bcfcomp);

    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    EnableIf_t<IsFabArray<MF>::value>
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

    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    EnableIf_t<IsFabArray<MF>::value>
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

#ifdef AMREX_USE_EB
    template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
    EnableIf_t<IsFabArray<MF>::value>
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

    template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
    EnableIf_t<IsFabArray<MF>::value>
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

    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    EnableIf_t<IsFabArray<MF>::value>
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

    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    EnableIf_t<IsFabArray<MF>::value>
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
}

#include <AMReX_FillPatchUtil_I.H>

#endif
