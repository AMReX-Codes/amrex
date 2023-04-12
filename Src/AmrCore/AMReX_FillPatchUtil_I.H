#ifndef AMREX_FillPatchUtil_I_H_
#define AMREX_FillPatchUtil_I_H_
#include <AMReX_Config.H>

namespace amrex {

namespace detail {

template <typename F, typename MF>
auto call_interp_hook (F const& f, MF& mf, int icomp, int ncomp)
    -> decltype(f(mf[0],Box(),icomp,ncomp))
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto& dfab = mf[mfi];
        const Box& dbx = dfab.box();
        f(dfab, dbx, icomp, ncomp);
    }
}

template <typename F, typename MF>
auto call_interp_hook (F const& f, MF& mf, int icomp, int ncomp)
    -> decltype(f(mf,icomp,ncomp))
{
    f(mf, icomp, ncomp);
}

}

template <typename Interp>
bool ProperlyNested (const IntVect& ratio, const IntVect& blocking_factor, int ngrow,
                     const IndexType& boxType, Interp* mapper)
{
    int ratio_max = ratio[0];
#if (AMREX_SPACEDIM > 1)
    ratio_max = std::max(ratio_max, ratio[1]);
#endif
#if (AMREX_SPACEDIM == 3)
    ratio_max = std::max(ratio_max, ratio[2]);
#endif
    // There are at least this many coarse cells outside fine grids
    // (except at physical boundaries).
    const IntVect& nbuf = blocking_factor / ratio_max;

    Box crse_box(IntVect(AMREX_D_DECL(0 ,0 ,0 )), IntVect(AMREX_D_DECL(4*nbuf[0]-1,
                                                                       4*nbuf[1]-1,
                                                                       4*nbuf[2]-1)));
    crse_box.convert(boxType);
    Box fine_box(nbuf, IntVect(AMREX_D_DECL(3*nbuf[0]-1,3*nbuf[1]-1,3*nbuf[2]-1)));
    fine_box.convert(boxType);
    fine_box.refine(ratio_max);
    fine_box.grow(ngrow);
    const Box& fine_box_coarsened = mapper->CoarseBox(fine_box, ratio_max);
    return crse_box.contains(fine_box_coarsened);
}

template <typename MF, typename BC>
std::enable_if_t<IsFabArray<MF>::value>
FillPatchSingleLevel (MF& mf, Real time,
                      const Vector<MF*>& smf, const Vector<Real>& stime,
                      int scomp, int dcomp, int ncomp,
                      const Geometry& geom,
                      BC& physbcf, int bcfcomp)
{
    FillPatchSingleLevel(mf, mf.nGrowVect(), time, smf, stime, scomp, dcomp, ncomp,
                         geom, physbcf, bcfcomp);
}

template <typename MF, typename BC>
std::enable_if_t<IsFabArray<MF>::value>
FillPatchSingleLevel (MF& mf, IntVect const& nghost, Real time,
                      const Vector<MF*>& smf, const Vector<Real>& stime,
                      int scomp, int dcomp, int ncomp,
                      const Geometry& geom,
                      BC& physbcf, int bcfcomp)
{
    BL_PROFILE("FillPatchSingleLevel");

    AMREX_ASSERT(scomp+ncomp <= smf[0]->nComp());
    AMREX_ASSERT(dcomp+ncomp <= mf.nComp());
    AMREX_ASSERT(smf.size() == stime.size());
    AMREX_ASSERT(!smf.empty());
    AMREX_ASSERT(nghost.allLE(mf.nGrowVect()));

    if (smf.size() == 1)
    {
        if (&mf == smf[0] && scomp == dcomp) {
            mf.FillBoundary(dcomp, ncomp, nghost, geom.periodicity());
        } else {
            mf.ParallelCopy(*smf[0], scomp, dcomp, ncomp, IntVect{0}, nghost, geom.periodicity());
        }
    }
    else if (smf.size() == 2)
    {
        BL_ASSERT(smf[0]->boxArray() == smf[1]->boxArray());
        MF raii;
        MF * dmf;
        int destcomp;
        bool sameba;
        if (mf.boxArray() == smf[0]->boxArray() &&
            mf.DistributionMap() == smf[0]->DistributionMap())
        {
            dmf = &mf;
            destcomp = dcomp;
            sameba = true;
        } else {
            raii.define(smf[0]->boxArray(), smf[0]->DistributionMap(), ncomp, 0,
                        MFInfo(), smf[0]->Factory());

            dmf = &raii;
            destcomp = 0;
            sameba = false;
        }

        if ((dmf != smf[0] && dmf != smf[1]) || scomp != dcomp)
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*dmf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Real t0 = stime[0];
                const Real t1 = stime[1];
                auto const sfab0 = smf[0]->array(mfi);
                auto const sfab1 = smf[1]->array(mfi);
                auto       dfab  = dmf->array(mfi);

                if (time == t0)
                {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
                    {
                        dfab(i,j,k,n+destcomp) = sfab0(i,j,k,n+scomp);
                    });
                }
                else if (time == t1)
                {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
                    {
                        dfab(i,j,k,n+destcomp) = sfab1(i,j,k,n+scomp);
                    });
                }
                else if (! amrex::almostEqual(t0,t1))
                {
                    Real alpha = (t1-time)/(t1-t0);
                    Real beta = (time-t0)/(t1-t0);
                    AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
                    {
                        dfab(i,j,k,n+destcomp) = alpha*sfab0(i,j,k,n+scomp)
                            +                     beta*sfab1(i,j,k,n+scomp);
                    });
                }
                else
                {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
                    {
                        dfab(i,j,k,n+destcomp) = sfab0(i,j,k,n+scomp);
                    });
                }
            }
        }

        if (sameba)
        {
            // Note that when sameba is true mf's BoxArray is nonoverlapping.
            // So FillBoundary is safe.
            mf.FillBoundary(dcomp, ncomp, nghost, geom.periodicity());
        }
        else
        {
            IntVect src_ngrow = IntVect::TheZeroVector();
            IntVect dst_ngrow = nghost;

            mf.ParallelCopy(*dmf, 0, dcomp, ncomp, src_ngrow, dst_ngrow, geom.periodicity());
        }
    }
    else {
        amrex::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
    }

    physbcf(mf, dcomp, ncomp, nghost, time, bcfcomp);
}

void FillPatchInterp (MultiFab& mf_fine_patch, int fcomp, MultiFab const& mf_crse_patch, int ccomp,
                      int ncomp, IntVect const& ng, const Geometry& cgeom, const Geometry& fgeom,
                      Box const& dest_domain, const IntVect& ratio,
                      MFInterpolater* mapper, const Vector<BCRec>& bcs, int bcscomp);

template <typename MF, typename Interp>
std::enable_if_t<IsFabArray<MF>::value && !std::is_same<Interp,MFInterpolater>::value>
FillPatchInterp (MF& mf_fine_patch, int fcomp, MF const& mf_crse_patch, int ccomp,
                 int ncomp, IntVect const& ng, const Geometry& cgeom, const Geometry& fgeom,
                 Box const& dest_domain, const IntVect& ratio,
                 Interp* mapper, const Vector<BCRec>& bcs, int bcscomp)
{
    BL_PROFILE("FillPatchInterp(Fab)");

    Box const& cdomain = amrex::convert(cgeom.Domain(), mf_fine_patch.ixType());
    int idummy=0;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        Vector<BCRec> bcr(ncomp);
        for (MFIter mfi(mf_fine_patch); mfi.isValid(); ++mfi)
        {
            auto& sfab = mf_crse_patch[mfi];
            const Box& sbx = sfab.box();

            auto& dfab = mf_fine_patch[mfi];
            Box const& dbx = amrex::grow(mfi.validbox(),ng) & dest_domain;

            amrex::setBC(sbx,cdomain,bcscomp,0,ncomp,bcs,bcr);
            mapper->interp(sfab, ccomp, dfab, fcomp, ncomp, dbx, ratio,
                           cgeom, fgeom, bcr, idummy, idummy, RunOn::Gpu);
        }
    }
}

template <typename MF>
std::enable_if_t<IsFabArray<MF>::value>
FillPatchInterp (MF& mf_fine_patch, int fcomp, MF const& mf_crse_patch, int ccomp,
                 int ncomp, IntVect const& ng, const Geometry& cgeom, const Geometry& fgeom,
                 Box const& dest_domain, const IntVect& ratio,
                 InterpBase* mapper, const Vector<BCRec>& bcs, int bcscomp)
{
    if (dynamic_cast<MFInterpolater*>(mapper)) {
        FillPatchInterp(mf_fine_patch, fcomp, mf_crse_patch, ccomp,
                        ncomp, ng, cgeom, fgeom, dest_domain, ratio,
                        static_cast<MFInterpolater*>(mapper), bcs, bcscomp);
    } else if (dynamic_cast<Interpolater*>(mapper)) {
        FillPatchInterp(mf_fine_patch, fcomp, mf_crse_patch, ccomp,
                        ncomp, ng, cgeom, fgeom, dest_domain, ratio,
                        static_cast<Interpolater*>(mapper), bcs, bcscomp);
    } else {
        amrex::Abort("FillPatchInterp: unknown InterpBase");
    }
}

template <typename MF, typename iMF, typename Interp>
std::enable_if_t<IsFabArray<MF>::value  && !std::is_same<Interp,MFInterpolater>::value>
InterpFace (Interp *interp,
            MF const& mf_crse_patch,    int crse_comp,
            MF&       mf_refined_patch, int fine_comp,
            int ncomp, const IntVect&   ratio,
            const iMF& solve_mask, const Geometry&  crse_geom, const Geometry&  fine_geom,
            int bcscomp, RunOn gpu_or_cpu,
            const Vector<BCRec>& bcs)
{
    Vector<BCRec> bcr(ncomp);
    Box const& cdomain = amrex::convert(crse_geom.Domain(), mf_crse_patch.ixType());
    for (MFIter mfi(mf_refined_patch);mfi.isValid(); ++mfi)
    {
        auto& sfab = mf_crse_patch[mfi];
        const Box& sbx = sfab.box();
        auto& dfab = mf_refined_patch[mfi];
        Box const& dbx = dfab.box();
        auto& ifab = solve_mask[mfi];
        amrex::setBC(sbx,cdomain,bcscomp,0,ncomp,bcs,bcr);
        interp->interp_face(sfab,crse_comp,dfab,fine_comp,ncomp,
                            dbx, ratio, ifab, crse_geom, fine_geom,
                            bcr, bcscomp, gpu_or_cpu);
    }
}

template <typename MF, typename iMF>
std::enable_if_t<IsFabArray<MF>::value>
InterpFace (InterpBase *interp,
            MF const& mf_crse_patch,    int crse_comp,
            MF&       mf_refined_patch, int fine_comp,
            int ncomp, const IntVect&   ratio,
            const iMF& solve_mask, const Geometry&  crse_geom, const Geometry&  fine_geom,
            int bccomp, RunOn gpu_or_cpu,
            const Vector<BCRec>& bcs)
{
    if (dynamic_cast<MFInterpolater*>(interp)){
        InterpFace(static_cast<MFInterpolater*>(interp),
                    mf_crse_patch, crse_comp,mf_refined_patch, fine_comp,
                    ncomp, ratio, solve_mask, crse_geom, fine_geom, bccomp,
                    gpu_or_cpu, bcs);
    }
    else if (dynamic_cast<Interpolater*>(interp)){
        InterpFace(static_cast<Interpolater*>(interp),
                    mf_crse_patch, crse_comp,mf_refined_patch, fine_comp,
                    ncomp, ratio, solve_mask, crse_geom, fine_geom, bccomp,
                    gpu_or_cpu, bcs);
    }
    else {
        amrex::Abort("InterpFace: unknown InterpBase");
    }
}


namespace {

// ======== FArrayBox

    template <typename MF,
              typename std::enable_if<std::is_same<typename MF::FABType::value_type,
                                                   FArrayBox>::value,
                                      int>::type = 0>
    MF make_mf_crse_patch (FabArrayBase::FPinfo const& fpc, int ncomp)
    {
        MF mf_crse_patch(fpc.ba_crse_patch, fpc.dm_patch, ncomp, 0, MFInfo(),
                         *fpc.fact_crse_patch);
        return mf_crse_patch;
    }

    template <typename MF,
              typename std::enable_if<std::is_same<typename MF::FABType::value_type,
                                                   FArrayBox>::value,
                                      int>::type = 0>
    MF make_mf_crse_patch (FabArrayBase::FPinfo const& fpc, int ncomp, IndexType idx_type)
    {
        MF mf_crse_patch(amrex::convert(fpc.ba_crse_patch, idx_type), fpc.dm_patch,
                         ncomp, 0, MFInfo(), *fpc.fact_crse_patch);
        return mf_crse_patch;
    }

    template <typename MF,
              typename std::enable_if<std::is_same<typename MF::FABType::value_type,
                                                   FArrayBox>::value,
                                      int>::type = 0>
    MF make_mf_fine_patch (FabArrayBase::FPinfo const& fpc, int ncomp)
    {
        MF mf_fine_patch(fpc.ba_fine_patch, fpc.dm_patch, ncomp, 0, MFInfo(),
                         *fpc.fact_fine_patch);
        return mf_fine_patch;
    }

    template <typename MF,
              typename std::enable_if<std::is_same<typename MF::FABType::value_type,
                                                   FArrayBox>::value,
                                      int>::type = 0>
    MF make_mf_fine_patch (FabArrayBase::FPinfo const& fpc, int ncomp, IndexType idx_type)
    {
        MF mf_fine_patch(amrex::convert(fpc.ba_fine_patch, idx_type), fpc.dm_patch,
                         ncomp, 0, MFInfo(), *fpc.fact_fine_patch);
        return mf_fine_patch;
    }

    template <typename MF,
              typename std::enable_if<std::is_same<typename MF::FABType::value_type,
                                                   FArrayBox>::value,
                                      int>::type = 0>
    MF make_mf_refined_patch (FabArrayBase::FPinfo const& fpc, int ncomp, IndexType idx_type, IntVect ratio)
    {
        MF mf_refined_patch(amrex::convert( amrex::refine( amrex::coarsen(fpc.ba_fine_patch, ratio), ratio), idx_type),
                            fpc.dm_patch, ncomp, 0, MFInfo(), *fpc.fact_fine_patch);
        return mf_refined_patch;
    }

    template <typename MF,
              typename std::enable_if<std::is_same<typename MF::FABType::value_type,
                                                   FArrayBox>::value,
                                      int>::type = 0>
    MF make_mf_crse_mask (FabArrayBase::FPinfo const& fpc, int ncomp, IndexType idx_type, IntVect ratio)
    {
        MF mf_crse_mask(amrex::convert(amrex::coarsen(fpc.ba_fine_patch, ratio), idx_type),
                        fpc.dm_patch, ncomp, 0, MFInfo(), *fpc.fact_fine_patch);
        return mf_crse_mask;
    }

    template <typename MF,
              typename std::enable_if<std::is_same<typename MF::FABType::value_type,
                                                   FArrayBox>::value,
                                      int>::type = 0>
    void mf_set_domain_bndry (MF &mf, Geometry const & geom)
    {
        mf.setDomainBndry(std::numeric_limits<Real>::quiet_NaN(), geom);
    }


// ======== Not FArrayBox

    template <typename MF,
              typename std::enable_if<!std::is_same<typename MF::FABType::value_type,
                                                    FArrayBox>::value,
                                      int>::type = 0>
    MF make_mf_crse_patch (FabArrayBase::FPinfo const& fpc, int ncomp)
    {
        return MF(fpc.ba_crse_patch, fpc.dm_patch, ncomp, 0);
    }

    template <typename MF,
              typename std::enable_if<!std::is_same<typename MF::FABType::value_type,
                                                    FArrayBox>::value,
                                      int>::type = 0>
    MF make_mf_crse_patch (FabArrayBase::FPinfo const& fpc, int ncomp, IndexType idx_type)
    {
        return MF(amrex::convert(fpc.ba_crse_patch, idx_type), fpc.dm_patch, ncomp, 0);
    }

    template <typename MF,
              typename std::enable_if<!std::is_same<typename MF::FABType::value_type,
                                                    FArrayBox>::value,
                                      int>::type = 0>
    MF make_mf_fine_patch (FabArrayBase::FPinfo const& fpc, int ncomp)
    {
        return MF(fpc.ba_fine_patch, fpc.dm_patch, ncomp, 0);
    }

    template <typename MF,
              typename std::enable_if<!std::is_same<typename MF::FABType::value_type,
                                                    FArrayBox>::value,
                                      int>::type = 0>
    MF make_mf_fine_patch (FabArrayBase::FPinfo const& fpc, int ncomp, IndexType idx_type)
    {
        return MF(amrex::convert(fpc.ba_fine_patch, idx_type), fpc.dm_patch, ncomp, 0);
    }

    template <typename MF,
              typename std::enable_if<!std::is_same<typename MF::FABType::value_type,
                                                   FArrayBox>::value,
                                      int>::type = 0>
    MF make_mf_refined_patch (FabArrayBase::FPinfo const& fpc, int ncomp, IndexType idx_type, IntVect ratio)
    {
        return MF(amrex::convert( amrex::refine( amrex::coarsen(fpc.ba_fine_patch, ratio), ratio), idx_type), fpc.dm_patch, ncomp, 0);
    }

    template <typename MF,
              typename std::enable_if<!std::is_same<typename MF::FABType::value_type,
                                                   FArrayBox>::value,
                                      int>::type = 0>
    MF make_mf_crse_mask (FabArrayBase::FPinfo const& fpc, int ncomp, IndexType idx_type, IntVect ratio)
    {
        return MF(amrex::convert(amrex::coarsen(fpc.ba_fine_patch, ratio), idx_type), fpc.dm_patch, ncomp, 0);
    }

    template <typename MF,
              typename std::enable_if<!std::is_same<typename MF::FABType::value_type,
                                                    FArrayBox>::value,
                                      int>::type = 0>
    void mf_set_domain_bndry (MF &/*mf*/, Geometry const & /*geom*/)
    {
        // nothing
    }

    template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels_doit (MF& mf, IntVect const& nghost, Real time,
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
                             const PostInterpHook& post_interp,
                             EB2::IndexSpace const* index_space)
    {
        BL_PROFILE("FillPatchTwoLevels");

        if (nghost.max() > 0 || mf.getBDKey() != fmf[0]->getBDKey())
        {
            const InterpolaterBoxCoarsener& coarsener = mapper->BoxCoarsener(ratio);

            // Test for Face-centered data
            if ( AMREX_D_TERM(  mf.ixType().nodeCentered(0),
                              + mf.ixType().nodeCentered(1),
                              + mf.ixType().nodeCentered(2) ) == 1 )
            {
                if ( !dynamic_cast<Interpolater*>(mapper) ){
                    amrex::Abort("This interpolater has not yet implemented a version for face-based data");
                }

                // Convert to cell-centered MF meta-data for FPInfo.
                MF mf_cc_dummy( amrex::convert(mf.boxArray(), IntVect::TheZeroVector()),
                                mf.DistributionMap(), ncomp, nghost, MFInfo().SetAlloc(false) );
                MF fmf_cc_dummy( amrex::convert(fmf[0]->boxArray(), IntVect::TheZeroVector()),
                                 fmf[0]->DistributionMap(), ncomp, nghost, MFInfo().SetAlloc(false) );

                const FabArrayBase::FPinfo& fpc = FabArrayBase::TheFPinfo(fmf_cc_dummy, mf_cc_dummy,
                                                                          nghost,
                                                                          coarsener,
                                                                          fgeom,
                                                                          cgeom,
                                                                          index_space);

                if ( ! fpc.ba_crse_patch.empty())
                {
                    MF mf_crse_patch     = make_mf_crse_patch<MF>      (fpc, ncomp, mf.boxArray().ixType());
                    // Must make sure fine exists under needed coarse faces.
                    // It stores values for the final (interior) interpolation,
                    // which is done from this fine MF that's been partially filled
                    // (with only faces overlying coarse having valid data).
                    MF mf_refined_patch  = make_mf_refined_patch<MF>   (fpc, ncomp, mf.boxArray().ixType(), ratio);
                    auto solve_mask = make_mf_crse_mask<iMultiFab>(fpc, ncomp, mf.boxArray().ixType(), ratio);

                    mf_set_domain_bndry(mf_crse_patch, cgeom);
                    FillPatchSingleLevel(mf_crse_patch, time, cmf, ct, scomp, 0, ncomp,
                                         cgeom, cbc, cbccomp);

                    mf_set_domain_bndry(mf_refined_patch, fgeom);
                    FillPatchSingleLevel(mf_refined_patch, time, fmf, ft, scomp, 0, ncomp,
                                         fgeom, fbc, fbccomp);

                    // Aliased MFs, used to allow CPC caching.
                    MF mf_known( amrex::coarsen(fmf[0]->boxArray(), ratio), fmf[0]->DistributionMap(),
                                 ncomp, nghost, MFInfo().SetAlloc(false) );
                    MF mf_solution( amrex::coarsen(mf_refined_patch.boxArray(), ratio), mf_refined_patch.DistributionMap(),
                                    ncomp, 0, MFInfo().SetAlloc(false) );

                    const FabArrayBase::CPC mask_cpc( mf_solution, IntVect::TheZeroVector(),
                                                      mf_known, IntVect::TheZeroVector(),
                                                      fgeom.periodicity());

                    solve_mask.setVal(1);                   // Values to solve.
                    solve_mask.setVal(0, mask_cpc, 0, 1);   // Known values.

                    detail::call_interp_hook(pre_interp, mf_crse_patch, 0, ncomp);

                    InterpFace(mapper, mf_crse_patch, 0, mf_refined_patch, 0, ncomp,
                               ratio, solve_mask, cgeom, fgeom, bcscomp, RunOn::Gpu, bcs);

                    detail::call_interp_hook(post_interp, mf_refined_patch, 0, ncomp);

                    bool aliasing = false;
                    for (auto const& fmf_a : fmf) {
                        aliasing = aliasing || (&mf == fmf_a);
                    }
                    if (aliasing) {
                        mf.ParallelCopyToGhost(mf_refined_patch, 0, dcomp, ncomp,
                                               IntVect{0}, nghost);
                    } else {
                        mf.ParallelCopy(mf_refined_patch, 0, dcomp, ncomp,
                                        IntVect{0}, nghost);
                    }
                }
            }
            else
            {
                const FabArrayBase::FPinfo& fpc = FabArrayBase::TheFPinfo(*fmf[0], mf,
                                                                          nghost,
                                                                          coarsener,
                                                                          fgeom,
                                                                          cgeom,
                                                                          index_space);

                if ( ! fpc.ba_crse_patch.empty())
                {

                    MF mf_crse_patch = make_mf_crse_patch<MF>(fpc, ncomp);
                    mf_set_domain_bndry (mf_crse_patch, cgeom);

                    FillPatchSingleLevel(mf_crse_patch, time, cmf, ct, scomp, 0, ncomp, cgeom, cbc, cbccomp);

                    MF mf_fine_patch = make_mf_fine_patch<MF>(fpc, ncomp);

                    detail::call_interp_hook(pre_interp, mf_crse_patch, 0, ncomp);

                    FillPatchInterp(mf_fine_patch, 0, mf_crse_patch, 0,
                                    ncomp, IntVect(0), cgeom, fgeom,
                                    amrex::grow(amrex::convert(fgeom.Domain(),mf.ixType()),nghost),
                                    ratio, mapper, bcs, bcscomp);

                    detail::call_interp_hook(post_interp, mf_fine_patch, 0, ncomp);

                    mf.ParallelCopy(mf_fine_patch, 0, dcomp, ncomp, IntVect{0}, nghost);
                }
            }
        }

        FillPatchSingleLevel(mf, nghost, time, fmf, ft, scomp, dcomp, ncomp,
                             fgeom, fbc, fbccomp);
    }

    template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels_doit (Array<MF*, AMREX_SPACEDIM> const& mf, IntVect const& nghost, Real time,
                             const Vector<Array<MF*, AMREX_SPACEDIM> >& cmf, const Vector<Real>& ct,
                             const Vector<Array<MF*, AMREX_SPACEDIM> >& fmf, const Vector<Real>& ft,
                             int scomp, int dcomp, int ncomp,
                             const Geometry& cgeom, const Geometry& fgeom,
                             Array<BC, AMREX_SPACEDIM>& cbc, const Array<int, AMREX_SPACEDIM>& cbccomp,
                             Array<BC, AMREX_SPACEDIM>& fbc, const Array<int, AMREX_SPACEDIM>& fbccomp,
                             const IntVect& ratio,
                             Interp* mapper,
                             const Array<Vector<BCRec>, AMREX_SPACEDIM>& bcs, const Array<int, AMREX_SPACEDIM>& bcscomp,
                             const PreInterpHook& pre_interp,
                             const PostInterpHook& post_interp,
                             EB2::IndexSpace const* index_space)
    {
        BL_PROFILE("FillPatchTwoLevels (Array<MF*>)");

        using FAB = typename MF::FABType::value_type;
        using iFAB = typename iMultiFab::FABType::value_type;

        AMREX_ASSERT(AMREX_D_TERM(mf[0]->ixType().nodeCentered(0),
                               && mf[1]->ixType().nodeCentered(1),
                               && mf[2]->ixType().nodeCentered(2)));

        // These need to be true: (ba[0] == ba[1] == ba[2]) & (dm[0] == dm[1] == dm[2]).
        // Debatable whether these are required, or will be enforced elsewhere prior to this func.
        AMREX_ASSERT(AMREX_D_TERM(true,
                               && BoxArray::SameRefs(mf[0]->boxArray(), mf[1]->boxArray()),
                               && BoxArray::SameRefs(mf[0]->boxArray(), mf[2]->boxArray())));
/*
        AMREX_ASSERT(AMREX_D_TERM(true,
                               && DistributionMapping::SameRefs(mf[0]->DistributionMap(), mf[1]->DistributionMap()),
                               && DistributionMapping::SameRefs(mf[0]->DistributionMap(), mf[2]->DistributionMap())));
*/


        // Test all of them?
        if (nghost.max() > 0 || mf[0]->getBDKey() != fmf[0][0]->getBDKey())
        {
            const InterpolaterBoxCoarsener& coarsener = mapper->BoxCoarsener(ratio);

            // Convert to cell-centered MF meta-data for FPInfo.
            MF mf_cc_dummy( amrex::convert(mf[0]->boxArray(), IntVect::TheZeroVector()),
                            mf[0]->DistributionMap(), ncomp, nghost, MFInfo().SetAlloc(false) );
            MF fmf_cc_dummy( amrex::convert(fmf[0][0]->boxArray(), IntVect::TheZeroVector()),
                             fmf[0][0]->DistributionMap(), ncomp, nghost, MFInfo().SetAlloc(false) );

            const FabArrayBase::FPinfo& fpc = FabArrayBase::TheFPinfo(fmf_cc_dummy, mf_cc_dummy,
                                                                      nghost,
                                                                      coarsener,
                                                                      fgeom,
                                                                      cgeom,
                                                                      index_space);

            if ( !fpc.ba_crse_patch.empty() )
            {
                Array<MF, AMREX_SPACEDIM> mf_crse_patch;
                Array<MF, AMREX_SPACEDIM> mf_refined_patch;
                Array<iMultiFab, AMREX_SPACEDIM> solve_mask;

                for (int d=0; d<AMREX_SPACEDIM; ++d)
                {
                    mf_crse_patch[d]    = make_mf_crse_patch<MF>      (fpc, ncomp, mf[d]->boxArray().ixType());
                    mf_refined_patch[d] = make_mf_refined_patch<MF>   (fpc, ncomp, mf[d]->boxArray().ixType(), ratio);
                    solve_mask[d]       = make_mf_crse_mask<iMultiFab>(fpc, ncomp, mf[d]->boxArray().ixType(), ratio);

                    mf_set_domain_bndry(mf_crse_patch[d], cgeom);
                    Vector<MF*> cmf_time;
                    for (const auto & mfab : cmf)
                        { cmf_time.push_back(mfab[d]); }

                    FillPatchSingleLevel(mf_crse_patch[d], time, cmf_time, ct, scomp, 0, ncomp,
                                         cgeom, cbc[d], cbccomp[d]);

                    mf_set_domain_bndry(mf_refined_patch[d], fgeom);
                    Vector<MF*> fmf_time;
                    for (const auto & mfab : fmf)
                        { fmf_time.push_back(mfab[d]); }

                    FillPatchSingleLevel(mf_refined_patch[d], time, fmf_time, ft, scomp, 0, ncomp,
                                         fgeom, fbc[d], fbccomp[d]);


                    // Aliased MFs, used to allow CPC caching.
                    MF mf_known( amrex::coarsen(fmf[0][d]->boxArray(), ratio), fmf[0][d]->DistributionMap(),
                                    ncomp, nghost, MFInfo().SetAlloc(false) );
                    MF mf_solution( amrex::coarsen(mf_refined_patch[d].boxArray(), ratio), mf_refined_patch[d].DistributionMap(),
                                    ncomp, 0, MFInfo().SetAlloc(false) );

                    const FabArrayBase::CPC mask_cpc( mf_solution, IntVect::TheZeroVector(),
                                                      mf_known, IntVect::TheZeroVector(),
                                                      fgeom.periodicity() );

                    solve_mask[d].setVal(1);                   // Values to solve.
                    solve_mask[d].setVal(0, mask_cpc, 0, 1);   // Known values.
                }

                int idummy=0;
#ifdef AMREX_USE_OMP
//              bool cc = fpc.ba_crse_patch.ixType().cellCentered();
                bool cc = false;        // can anything be done to allow threading, or can the OpenMP just be removed?
#pragma omp parallel if (cc && Gpu::notInLaunchRegion() )
#endif
                {
                    Vector<Array<BCRec, AMREX_SPACEDIM> > bcr(ncomp);
                    for (MFIter mfi(mf_refined_patch[0]); mfi.isValid(); ++mfi)
                    {
                        Array<FAB*, AMREX_SPACEDIM> sfab{ AMREX_D_DECL( &(mf_crse_patch[0][mfi]),
                                                                        &(mf_crse_patch[1][mfi]),
                                                                        &(mf_crse_patch[2][mfi])  )};
                        Array<FAB*, AMREX_SPACEDIM> dfab{ AMREX_D_DECL( &(mf_refined_patch[0][mfi]),
                                                                        &(mf_refined_patch[1][mfi]),
                                                                        &(mf_refined_patch[2][mfi])  )};
                        Array<iFAB*, AMREX_SPACEDIM> mfab{ AMREX_D_DECL( &(solve_mask[0][mfi]),
                                                                         &(solve_mask[1][mfi]),
                                                                         &(solve_mask[2][mfi])  )};

                        const Box& sbx_cc = amrex::convert(sfab[0]->box(), IntVect::TheZeroVector());
                        const Box& dbx_cc = amrex::convert(dfab[0]->box(), IntVect::TheZeroVector());

                        for (int d=0; d<AMREX_SPACEDIM; ++d)
                        {
                            Vector<BCRec> bcr_d(ncomp);

                            amrex::setBC(sfab[d]->box(), amrex::convert(cgeom.Domain(), mf[d]->ixType()),
                                         bcscomp[d],0,ncomp,bcs[d],bcr_d);

                            for (int n=0; n<ncomp; ++n)
                                { bcr[n][d] = bcr_d[n]; }
                        }

                        pre_interp(sfab, sbx_cc, 0, ncomp);

                        mapper->interp_arr(sfab, 0, dfab, 0, ncomp, dbx_cc, ratio, mfab,
                                           cgeom, fgeom, bcr, idummy, idummy, RunOn::Gpu);

                        post_interp(dfab, dbx_cc, 0, ncomp);
                    }
                }

                for (int d=0; d<AMREX_SPACEDIM; ++d)
                {
                    bool aliasing = false;
                    for (auto const& fmf_a : fmf) {
                        aliasing = aliasing || (mf[d] == fmf_a[d]);
                    }
                    if (aliasing) {
                        mf[d]->ParallelCopyToGhost(mf_refined_patch[d], 0, dcomp, ncomp,
                                                   IntVect{0}, nghost);
                    } else {
                        mf[d]->ParallelCopy(mf_refined_patch[d], 0, dcomp, ncomp,
                                            IntVect{0}, nghost);
                    }
                }
            }
        }

        for (int d=0; d<AMREX_SPACEDIM; ++d)
        {
            Vector<MF*> fmf_time;
            for (auto const& ffab : fmf)
                 { fmf_time.push_back(ffab[d]); }

            FillPatchSingleLevel(*mf[d], nghost, time, fmf_time, ft, scomp, dcomp, ncomp,
                                 fgeom, fbc[d], fbccomp[d]);
        }
    }

} // Anonymous namespace

template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
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
                    const PreInterpHook& pre_interp,
                    const PostInterpHook& post_interp)
{
#ifdef AMREX_USE_EB
    EB2::IndexSpace const* index_space = EB2::TopIndexSpaceIfPresent();
#else
    EB2::IndexSpace const* index_space = nullptr;
#endif
    FillPatchTwoLevels_doit(mf,nghost,time,cmf,ct,fmf,ft,
                            scomp,dcomp,ncomp,cgeom,fgeom,
                            cbc,cbccomp,fbc,fbccomp,ratio,mapper,bcs,bcscomp,
                            pre_interp,post_interp,index_space);
}

template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
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
                    const PreInterpHook& pre_interp,
                    const PostInterpHook& post_interp)
{
#ifdef AMREX_USE_EB
    EB2::IndexSpace const* index_space = EB2::TopIndexSpaceIfPresent();
#else
    EB2::IndexSpace const* index_space = nullptr;
#endif

    FillPatchTwoLevels_doit(mf,mf.nGrowVect(),time,cmf,ct,fmf,ft,
                            scomp,dcomp,ncomp,cgeom,fgeom,
                            cbc,cbccomp,fbc,fbccomp,ratio,mapper,bcs,bcscomp,
                            pre_interp,post_interp,index_space);
}

template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
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
                    const PreInterpHook& pre_interp,
                    const PostInterpHook& post_interp)
{
#ifdef AMREX_USE_EB
    EB2::IndexSpace const* index_space = EB2::TopIndexSpaceIfPresent();
#else
    EB2::IndexSpace const* index_space = nullptr;
#endif

    FillPatchTwoLevels_doit(mf,nghost,time,cmf,ct,fmf,ft,
                            scomp,dcomp,ncomp,cgeom,fgeom,
                            cbc,cbccomp,fbc,fbccomp,ratio,mapper,bcs,bcscomp,
                            pre_interp,post_interp,index_space);
}

template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
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
                    const PreInterpHook& pre_interp,
                    const PostInterpHook& post_interp)
{
#ifdef AMREX_USE_EB
    EB2::IndexSpace const* index_space = EB2::TopIndexSpaceIfPresent();
#else
    EB2::IndexSpace const* index_space = nullptr;
#endif

    Array<int, AMREX_SPACEDIM> cbccomp_arr = {AMREX_D_DECL(cbccomp,cbccomp,cbccomp)};
    Array<int, AMREX_SPACEDIM> fbccomp_arr = {AMREX_D_DECL(fbccomp,fbccomp,fbccomp)};
    Array<int, AMREX_SPACEDIM> bcscomp_arr = {AMREX_D_DECL(bcscomp,bcscomp,bcscomp)};

    FillPatchTwoLevels_doit(mf,nghost,time,cmf,ct,fmf,ft,
                            scomp,dcomp,ncomp,cgeom,fgeom,
                            cbc,cbccomp_arr,fbc,fbccomp_arr,ratio,mapper,bcs,bcscomp_arr,
                            pre_interp,post_interp,index_space);
}

template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
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
                    const PreInterpHook& pre_interp,
                    const PostInterpHook& post_interp)
{
#ifdef AMREX_USE_EB
    EB2::IndexSpace const* index_space = EB2::TopIndexSpaceIfPresent();
#else
    EB2::IndexSpace const* index_space = nullptr;
#endif

    Array<int, AMREX_SPACEDIM> cbccomp_arr = {AMREX_D_DECL(cbccomp,cbccomp,cbccomp)};
    Array<int, AMREX_SPACEDIM> fbccomp_arr = {AMREX_D_DECL(fbccomp,fbccomp,fbccomp)};
    Array<int, AMREX_SPACEDIM> bcscomp_arr = {AMREX_D_DECL(bcscomp,bcscomp,bcscomp)};

    FillPatchTwoLevels_doit(mf,mf[0]->nGrowVect(),time,cmf,ct,fmf,ft,
                            scomp,dcomp,ncomp,cgeom,fgeom,
                            cbc,cbccomp_arr,fbc,fbccomp_arr,ratio,mapper,bcs,bcscomp_arr,
                            pre_interp,post_interp,index_space);
}

#ifdef AMREX_USE_EB
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
                    const PostInterpHook& post_interp)
{
    FillPatchTwoLevels_doit(mf,nghost,time,cmf,ct,fmf,ft,
                            scomp,dcomp,ncomp,cgeom,fgeom,
                            cbc,cbccomp,fbc,fbccomp,ratio,mapper,bcs,bcscomp,
                            pre_interp,post_interp,&index_space);
}

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
                    const PostInterpHook& post_interp)
{
    FillPatchTwoLevels_doit(mf,mf.nGrowVect(),time,cmf,ct,fmf,ft,
                            scomp,dcomp,ncomp,cgeom,fgeom,
                            cbc,cbccomp,fbc,fbccomp,ratio,mapper,bcs,bcscomp,
                            pre_interp,post_interp,&index_space);
}
#endif

template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
std::enable_if_t<IsFabArray<MF>::value>
InterpFromCoarseLevel (MF& mf, Real time,
                       const MF& cmf, int scomp, int dcomp, int ncomp,
                       const Geometry& cgeom, const Geometry& fgeom,
                       BC& cbc, int cbccomp,
                       BC& fbc, int fbccomp,
                       const IntVect& ratio,
                       Interp* mapper,
                       const Vector<BCRec>& bcs, int bcscomp,
                       const PreInterpHook& pre_interp,
                       const PostInterpHook& post_interp)
{
    InterpFromCoarseLevel(mf,mf.nGrowVect(),time,cmf,scomp,dcomp,ncomp,cgeom,fgeom,
                          cbc,cbccomp,fbc,fbccomp,ratio,mapper,bcs,bcscomp,
                          pre_interp,post_interp);
}

template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
std::enable_if_t<IsFabArray<MF>::value>
InterpFromCoarseLevel (Array<MF*, AMREX_SPACEDIM> const& mf, Real time,
                       const Array<MF*, AMREX_SPACEDIM>& cmf, int scomp, int dcomp, int ncomp,
                       const Geometry& cgeom, const Geometry& fgeom,
                       Array<BC, AMREX_SPACEDIM>& cbc, int cbccomp,
                       Array<BC, AMREX_SPACEDIM>& fbc, int fbccomp,
                       const IntVect& ratio,
                       Interp* mapper,
                       const Array<Vector<BCRec>, AMREX_SPACEDIM>& bcs, int bcscomp,
                       const PreInterpHook& pre_interp,
                       const PostInterpHook& post_interp)
{
    InterpFromCoarseLevel(mf,mf[0]->nGrowVect(),time,cmf,scomp,dcomp,ncomp,cgeom,fgeom,
                          cbc,cbccomp,fbc,fbccomp,ratio,mapper,bcs,bcscomp,
                          pre_interp,post_interp);
}

template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
std::enable_if_t<IsFabArray<MF>::value>
InterpFromCoarseLevel (MF& mf, IntVect const& nghost, Real time,
                       const MF& cmf, int scomp, int dcomp, int ncomp,
                       const Geometry& cgeom, const Geometry& fgeom,
                       BC& cbc, int cbccomp,
                       BC& fbc, int fbccomp,
                       const IntVect& ratio,
                       Interp* mapper,
                       const Vector<BCRec>& bcs, int bcscomp,
                       const PreInterpHook& pre_interp,
                       const PostInterpHook& post_interp)
{
    using FAB = typename MF::FABType::value_type;

    const InterpolaterBoxCoarsener& coarsener = mapper->BoxCoarsener(ratio);

    const BoxArray& ba = mf.boxArray();
    const DistributionMapping& dm = mf.DistributionMap();

    const IndexType& typ = ba.ixType();

    BL_ASSERT(typ == cmf.boxArray().ixType());

    Box fdomain_g( amrex::convert(fgeom.Domain(),mf.ixType()) );
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (fgeom.isPeriodic(i)) {
            fdomain_g.grow(i,nghost[i]);
        }
    }

    BoxArray ba_crse_patch(ba.size());
    {  // TODO: later we might want to cache this
        for (int i = 0, N = ba.size(); i < N; ++i)
        {
            Box bx = amrex::convert(amrex::grow(ba[i],nghost), typ);
            bx &= fdomain_g;
            ba_crse_patch.set(i, coarsener.doit(bx));
        }
    }

    MF mf_crse_patch;
#ifdef AMREX_USE_EB
    if (EB2::TopIndexSpaceIfPresent()) {
        auto factory = makeEBFabFactory(cgeom, ba_crse_patch, dm, {0,0,0}, EBSupport::basic);
        mf_crse_patch.define(ba_crse_patch, dm, ncomp, 0, MFInfo(), *factory);
    } else
#endif
    {
        mf_crse_patch.define(ba_crse_patch, dm, ncomp, 0);
    }
    mf_set_domain_bndry (mf_crse_patch, cgeom);

    mf_crse_patch.ParallelCopy(cmf, scomp, 0, ncomp, cgeom.periodicity());

    cbc(mf_crse_patch, 0, ncomp, mf_crse_patch.nGrowVect(), time, cbccomp);

    detail::call_interp_hook(pre_interp, mf_crse_patch, 0, ncomp);

    FillPatchInterp(mf, dcomp, mf_crse_patch, 0, ncomp, nghost, cgeom, fgeom, fdomain_g,
                    ratio, mapper, bcs, bcscomp);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FAB& dfab   = mf[mfi];
        Box dfab_bx = dfab.box();
        dfab_bx.grow(nghost-mf.nGrowVect());
        const Box& dbx = dfab_bx & fdomain_g;

        post_interp(dfab, dbx, dcomp, ncomp);
    }

    fbc(mf, dcomp, ncomp, nghost, time, fbccomp);
}

template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
std::enable_if_t<IsFabArray<MF>::value>
InterpFromCoarseLevel (Array<MF*, AMREX_SPACEDIM> const& mf, IntVect const& nghost, Real time,
                       const Array<MF*, AMREX_SPACEDIM>& cmf, int scomp, int dcomp, int ncomp,
                       const Geometry& cgeom, const Geometry& fgeom,
                       Array<BC, AMREX_SPACEDIM>& cbc, int cbccomp,
                       Array<BC, AMREX_SPACEDIM>& fbc, int fbccomp,
                       const IntVect& ratio,
                       Interp* mapper,
                       const Array<Vector<BCRec>, AMREX_SPACEDIM>& bcs, int bcscomp,
                       const PreInterpHook& pre_interp,
                       const PostInterpHook& post_interp)
{
    using FAB = typename MF::FABType::value_type;
    using iFAB = typename iMultiFab::FABType::value_type;

    const InterpolaterBoxCoarsener& coarsener = mapper->BoxCoarsener(ratio);
    const BoxArray& ba = mf[0]->boxArray();
    const DistributionMapping& dm = mf[0]->DistributionMap();

    AMREX_ASSERT(AMREX_D_TERM(mf[0]->ixType().nodeCentered(0),
                           && mf[1]->ixType().nodeCentered(1),
                           && mf[2]->ixType().nodeCentered(2)));

    // These need to be true: (ba[0] == ba[1] == ba[2]) & (dm[0] == dm[1] == dm[2]).
    // Debatable whether these are required, or will be enforced elsewhere prior to this func.
    AMREX_ASSERT(AMREX_D_TERM(true,
                           && BoxArray::SameRefs(mf[0]->boxArray(), mf[1]->boxArray()),
                           && BoxArray::SameRefs(mf[0]->boxArray(), mf[2]->boxArray())));
/*
    AMREX_ASSERT(AMREX_D_TERM(true,
                           && DistributionMapping::SameRefs(mf[0]->DistributionMap(), mf[1]->DistributionMap()),
                           && DistributionMapping::SameRefs(mf[0]->DistributionMap(), mf[2]->DistributionMap())));
*/

    // If needed, adjust to fully overlap the coarse cells.
    IntVect nghost_adj = nghost;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        if (nghost[d] % ratio[d] != 0) {
            nghost_adj[d] += ratio[d] - (nghost[d] % ratio[d]);
        }
    }

    Array<MF*, AMREX_SPACEDIM> mf_local = mf;
    int dcomp_adj = dcomp;
    Array<std::unique_ptr<MF>, AMREX_SPACEDIM> mf_temp;
    if (! nghost.allGE(nghost_adj)) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            mf_temp[d] = std::make_unique<MF>(mf[d]->boxArray(),
                                              mf[d]->DistributionMap(), ncomp, nghost_adj);
            mf_local[d] = mf_temp[d].get();
        }
        dcomp_adj = 0;
    }

    // Create a cell-centered boxArray of the region to interp.
    // Convert this boxArray and domain as needed.
    Box fdomain = amrex::convert(fgeom.Domain(), IntVect::TheZeroVector());
    Box fdomain_g(fdomain);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (fgeom.isPeriodic(d)) {
            fdomain_g.grow(d,nghost_adj[d]);
        }
    }

    // Build patches, using domain to account for periodic bcs.
    BoxArray ba_crse_patch(ba.size());
    {  // TODO: later we might want to cache this
        for (int i = 0, N = ba.size(); i < N; ++i)
        {
            Box bx = amrex::convert(amrex::grow(ba[i], nghost_adj), IntVect::TheZeroVector());
            bx &= fdomain_g;
            ba_crse_patch.set(i, coarsener.doit(bx));
        }
    }

    Array<MF, AMREX_SPACEDIM> mf_crse_patch;
    for (int d = 0; d<AMREX_SPACEDIM; ++d)
    {
        IndexType typ = mf[d]->boxArray().ixType();
        BoxArray ba_crse_idxed = amrex::convert(ba_crse_patch, typ);

#ifdef AMREX_USE_EB
        auto crse_factory = makeEBFabFactory(cgeom, ba_crse_idxed, dm, {0,0,0}, EBSupport::basic);
        mf_crse_patch[d].define(ba_crse_idxed, dm, ncomp, 0, MFInfo(), *crse_factory);
#else
        mf_crse_patch[d].define(ba_crse_idxed, dm, ncomp, 0);
#endif
        mf_set_domain_bndry(mf_crse_patch[d], cgeom);

        mf_crse_patch[d].ParallelCopy(*(cmf[d]), scomp, 0, ncomp, cgeom.periodicity());
        cbc[d](mf_crse_patch[d], 0, ncomp, mf_crse_patch[d].nGrowVect(), time, cbccomp);
    }

    int idummy1=0, idummy2=0;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        Vector<Array<BCRec, AMREX_SPACEDIM> > bcr(ncomp);

        // Empty containers describing that all points must be solved (no mask).
        Array<iFAB*, AMREX_SPACEDIM> mfab{ AMREX_D_DECL( nullptr, nullptr, nullptr ) };

        for (MFIter mfi(mf_crse_patch[0]); mfi.isValid(); ++mfi)
        {
            Array<FAB*, AMREX_SPACEDIM>  sfab{ AMREX_D_DECL( &(mf_crse_patch[0][mfi]),
                                                             &(mf_crse_patch[1][mfi]),
                                                             &(mf_crse_patch[2][mfi]) )};
            Array<FAB*, AMREX_SPACEDIM>  dfab{ AMREX_D_DECL( &(*mf_local[0])[mfi],
                                                             &(*mf_local[1])[mfi],
                                                             &(*mf_local[2])[mfi] )};

            const Box& sbx_cc = amrex::convert(sfab[0]->box(), IntVect::TheZeroVector());
            Box dfab_cc = amrex::convert(dfab[0]->box(), IntVect::TheZeroVector());
            const Box& dbx_cc = dfab_cc & fdomain_g;

            for (int d=0; d<AMREX_SPACEDIM; ++d)
            {
                Vector<BCRec> bcr_d(ncomp);

                amrex::setBC(sfab[d]->box(),
                             amrex::convert(cgeom.Domain(), sfab[d]->box().ixType()),
                             bcscomp,0,ncomp,bcs[d],bcr_d);

                for (int n=0; n<ncomp; ++n)
                    { bcr[n][d] = bcr_d[n]; }
            }

            pre_interp(sfab, sbx_cc, 0, ncomp);

            mapper->interp_arr(sfab, 0, dfab, 0, ncomp, dbx_cc, ratio, mfab,
                               cgeom, fgeom, bcr, idummy1, idummy2, RunOn::Gpu);

            post_interp(dfab, dbx_cc, 0, ncomp);
        }
    }

    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
        if (mf[d] != mf_local[d]) {
            amrex::Copy(*mf[d], *mf_local[d], 0, dcomp_adj, ncomp, nghost);
        }

        fbc[d](*mf[d], dcomp, ncomp, nghost, time, fbccomp);
    }
}

}

#endif
