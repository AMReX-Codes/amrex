#ifndef AMREX_FillPatchUtil_H_
#define AMREX_FillPatchUtil_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Interpolater.H>
#include <array>

namespace amrex
{
    class InterpHook {
    public:
        virtual ~InterpHook () {}
        virtual void operator() (FArrayBox& fab, const Box& bx, int icomp, int ncomp) const = 0;
    };

    class NullInterpHook final
        : public InterpHook
    {
    public:
        virtual void operator() (FArrayBox& fab, const Box& bx, int icomp, int ncomp) const final {}
    };

    bool ProperlyNested (const IntVect& ratio, const IntVect& blockint_factor, int ngrow,
			 const IndexType& boxType, Interpolater* mapper);

    void FillPatchSingleLevel (MultiFab& mf, Real time,
			       const Vector<MultiFab*>& smf, const Vector<Real>& stime,
			       int scomp, int dcomp, int ncomp,
			       const Geometry& geom, PhysBCFunctBase& physbcf);

    // Specialization for non-interpolating in time FillPatch. NOTE: the time
    // variable is used by the PhysBC
    void FillPatchSingleLevel (MultiFab & mf, Real time, const MultiFab & smf,
                               int scomp, int dcomp, int ncomp,
                               const Geometry & geom, PhysBCFunctBase & physbcf);

    void FillPatchTwoLevels (MultiFab& mf, Real time,
			     const Vector<MultiFab*>& cmf, const Vector<Real>& ct,
			     const Vector<MultiFab*>& fmf, const Vector<Real>& ft,
			     int scomp, int dcomp, int ncomp,
			     const Geometry& cgeom, const Geometry& fgeom,
			     PhysBCFunctBase& cbc, PhysBCFunctBase& fbc,
			     const IntVect& ratio,
			     Interpolater* mapper, const BCRec& bcs,
                             const InterpHook& pre_interp = NullInterpHook(),
                             const InterpHook& post_interp = NullInterpHook());

    void FillPatchTwoLevels (MultiFab& mf, Real time,
			     const Vector<MultiFab*>& cmf, const Vector<Real>& ct,
			     const Vector<MultiFab*>& fmf, const Vector<Real>& ft,
			     int scomp, int dcomp, int ncomp,
			     const Geometry& cgeom, const Geometry& fgeom,
			     PhysBCFunctBase& cbc, PhysBCFunctBase& fbc,
			     const IntVect& ratio,
			     Interpolater* mapper, const Vector<BCRec>& bcs,
                             const InterpHook& pre_interp = NullInterpHook(),
                             const InterpHook& post_interp = NullInterpHook());

    // Specialization for non-interpolating in time FillPatch. NOTE: the time
    // variable is used by the PhysBC
    void FillPatchTwoLevels (MultiFab& mf, Real time,
                             const MultiFab & cmf, const MultiFab & fmf,
                             int scomp, int dcomp, int ncomp,
                             const Geometry& cgeom, const Geometry& fgeom,
                             PhysBCFunctBase& cbc, PhysBCFunctBase& fbc,
                             const IntVect& ratio,
                             Interpolater* mapper, const BCRec& bcs,
                             const InterpHook& pre_interp = NullInterpHook(),
                             const InterpHook& post_interp = NullInterpHook());

    // Specialization for non-interpolating in time FillPatch. NOTE: the time
    // variable is used by the PhysBC
    void FillPatchTwoLevels (MultiFab& mf, Real time,
                             const MultiFab & cmf, const MultiFab & fmf,
                             int scomp, int dcomp, int ncomp,
                             const Geometry& cgeom, const Geometry& fgeom,
                             PhysBCFunctBase& cbc, PhysBCFunctBase& fbc,
                             const IntVect& ratio,
                             Interpolater* mapper, const Vector<BCRec>& bcs,
                             const InterpHook& pre_interp = NullInterpHook(),
                             const InterpHook& post_interp = NullInterpHook());

    void InterpFromCoarseLevel (MultiFab& mf, Real time,
				const MultiFab& cmf, int scomp, int dcomp, int ncomp,
				const Geometry& cgeom, const Geometry& fgeom,
				PhysBCFunctBase& cbc, PhysBCFunctBase& fbc, const IntVect& ratio,
				Interpolater* mapper, const BCRec& bcs,
                                const InterpHook& pre_interp = NullInterpHook(),
                                const InterpHook& post_interp = NullInterpHook());


    void InterpFromCoarseLevel (MultiFab& mf, Real time,
				const MultiFab& cmf, int scomp, int dcomp, int ncomp,
				const Geometry& cgeom, const Geometry& fgeom,
				PhysBCFunctBase& cbc, PhysBCFunctBase& fbc, const IntVect& ratio,
				Interpolater* mapper, const Vector<BCRec>& bcs,
                                const InterpHook& pre_interp = NullInterpHook(),
                                const InterpHook& post_interp = NullInterpHook());

    void InterpCrseFineBndryBfield (const std::array<MultiFab,AMREX_SPACEDIM>& crse,
                                    std::array<MultiFab,AMREX_SPACEDIM>& fine,
                                    const Geometry& cgeom, const Geometry& fgeom,
                                    int ref_ratio);

    enum InterpEM_t { InterpE, InterpB};

    void InterpCrseFineBndryEMfield (InterpEM_t interp_type,
                                     const std::array<MultiFab,AMREX_SPACEDIM>& crse,
                                     std::array<MultiFab,AMREX_SPACEDIM>& fine,
                                     const Geometry& cgeom, const Geometry& fgeom,
                                     int ref_ratio);
}

#endif
