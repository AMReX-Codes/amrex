#ifndef AMREX_FillPatchUtil_I_H_
#define AMREX_FillPatchUtil_I_H_

namespace amrex {

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
EnableIf_t<IsFabArray<MF>::value>
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
EnableIf_t<IsFabArray<MF>::value>
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
    AMREX_ASSERT(smf.size() != 0);
    AMREX_ASSERT(nghost.allLE(mf.nGrowVect()));

    if (smf.size() == 1)
    {
        if (&mf == smf[0] and scomp == dcomp) {
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
        if (mf.boxArray() == smf[0]->boxArray() and
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

        if ((dmf != smf[0] and dmf != smf[1]) or scomp != dcomp)
        {
#ifdef _OPENMP
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
                else if (std::abs(t1-t0) > 1.e-16)
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

namespace {
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
              typename std::enable_if<std::is_same<typename MF::FABType::value_type,
                                                   FArrayBox>::value,
                                      int>::type = 0>
    void mf_set_domain_bndry (MF &mf, Geometry const & geom)
    {
        mf.setDomainBndry(std::numeric_limits<Real>::quiet_NaN(), geom);
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
    EnableIf_t<IsFabArray<MF>::value>
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

        using FAB = typename MF::FABType::value_type;

        if (nghost.max() > 0 || mf.getBDKey() != fmf[0]->getBDKey())
        {
            const InterpolaterBoxCoarsener& coarsener = mapper->BoxCoarsener(ratio);

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

                Box const& fdomain = amrex::convert(fgeom.Domain(),mf.ixType());
                int idummy=0;
#ifdef _OPENMP
                bool cc = fpc.ba_crse_patch.ixType().cellCentered();
#pragma omp parallel if (cc && Gpu::notInLaunchRegion())
#endif
                {
                    Vector<BCRec> bcr(ncomp);
                    for (MFIter mfi(mf_fine_patch); mfi.isValid(); ++mfi)
                    {
                        FAB& sfab = mf_crse_patch[mfi];
                        FAB& dfab = mf_fine_patch[mfi];
                        const Box& dbx = dfab.box();

                        amrex::setBC(dbx,fdomain,bcscomp,0,ncomp,bcs,bcr);

                        pre_interp(sfab, sfab.box(), 0, ncomp);

                        mapper->interp(sfab, 0, dfab, 0, ncomp, dbx, ratio,
                                       cgeom, fgeom, bcr, dcomp, idummy, RunOn::Gpu);

                        post_interp(dfab, dbx, 0, ncomp);
                    }
                }

                mf.ParallelCopy(mf_fine_patch, 0, dcomp, ncomp, IntVect{0}, nghost);
	    }
	}

        FillPatchSingleLevel(mf, nghost, time, fmf, ft, scomp, dcomp, ncomp,
                             fgeom, fbc, fbccomp);
    }

    template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
    EnableIf_t<IsFabArray<MF>::value>
    FillPatchTwoLevels_doit (Array<MF*, AMREX_SPACEDIM> const& mf, IntVect const& nghost, Real time,
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
                             const PostInterpHook& post_interp,
                             EB2::IndexSpace const* index_space)
    {
	BL_PROFILE("FillPatchTwoLevels-SPACEDIM");

        // Appropriate ASSERTs here?
        // BoxArray's equal except for IndexType?
        // IndexType?

        using FAB = typename MF::FABType::value_type;

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

            amrex::Print() << "fpc.ba_c = " << fpc.ba_crse_patch << std::endl;
            amrex::Print() << "fpc.ba_f = " << fpc.ba_fine_patch << std::endl;

	    if ( !fpc.ba_crse_patch.empty() )
	    {
                Array<MF, AMREX_SPACEDIM> mf_crse_patch;
                Array<MF, AMREX_SPACEDIM> mf_fine_patch;
                for (int i=0; i<AMREX_SPACEDIM; ++i)
                {
                    mf_crse_patch[i] = make_mf_crse_patch<MF>(fpc, ncomp, mf[i]->boxArray().ixType());
                    mf_fine_patch[i] = make_mf_fine_patch<MF>(fpc, ncomp, mf[i]->boxArray().ixType());
                }

                for (int i=0; i<AMREX_SPACEDIM; ++i)
                {
                    mf_set_domain_bndry(mf_crse_patch[i], cgeom);

                    Vector<MF*> cmf_time;
                    for (const auto & cfab : cmf)
                        { cmf_time.push_back(cfab[i]); }

                    FillPatchSingleLevel(mf_crse_patch[i], time, cmf_time, ct, scomp, 0, ncomp,
                                         cgeom, cbc[i], cbccomp);
                }

                int idummy=0;
#ifdef _OPENMP
                bool cc = fpc.ba_crse_patch.ixType().cellCentered();   // What do we do here?
#pragma omp parallel if (cc && Gpu::notInLaunchRegion() )
#endif
                {
                    Vector<Array<BCRec, AMREX_SPACEDIM> > bcr(ncomp);
                    for (MFIter mfi(mf_fine_patch[0]); mfi.isValid(); ++mfi)
                    {
                        Array<FAB*, AMREX_SPACEDIM> sfab{ AMREX_D_DECL( &(mf_crse_patch[0][mfi]),
                                                                        &(mf_crse_patch[1][mfi]),
                                                                        &(mf_crse_patch[2][mfi]) )};
                        Array<FAB*, AMREX_SPACEDIM> dfab{ AMREX_D_DECL( &(mf_fine_patch[0][mfi]),
                                                                        &(mf_fine_patch[1][mfi]),
                                                                        &(mf_fine_patch[2][mfi]) )};

                        const Box& sbx_cc = amrex::convert(sfab[0]->box(), IntVect::TheZeroVector());
                        const Box& dbx_cc = amrex::convert(dfab[0]->box(), IntVect::TheZeroVector());

                        for (int d=0; d<AMREX_SPACEDIM; ++d)
                        {
                            Vector<BCRec> bcr_d(ncomp);

                            amrex::setBC(dfab[d]->box(), amrex::convert(fgeom.Domain(), mf[d]->boxArray().ixType()),
                                         bcscomp,0,ncomp,bcs[d],bcr_d);

                            for (int n=0; n<ncomp; ++n)
                                { bcr[n][d] = bcr_d[n]; }
                        }

                        pre_interp(sfab, sbx_cc, 0, ncomp);

                        mapper->interp_arr(sfab, 0, dfab, 0, ncomp, dbx_cc, ratio,
                                           cgeom, fgeom, bcr, dcomp, idummy, RunOn::Gpu);

                        post_interp(dfab, dbx_cc, 0, ncomp);
                    }
                }

                for (int i=0; i<AMREX_SPACEDIM; ++i)
                {
                    mf[i]->ParallelCopy(mf_fine_patch[i], 0, dcomp, ncomp, IntVect{0}, nghost);
                }
	    }
	}

        for (int i=0; i<AMREX_SPACEDIM; ++i)
        {
            Vector<MF*> fmf_time;
            for (auto ffab : fmf)
                 { fmf_time.push_back(ffab[i]); }

	    FillPatchSingleLevel(*mf[i], nghost, time, fmf_time, ft, scomp, dcomp, ncomp,
                                 fgeom, fbc[i], fbccomp);
        }
    }

} // Anonymous namespace

template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
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
EnableIf_t<IsFabArray<MF>::value>
FillPatchTwoLevels (Array<MF*, AMREX_SPACEDIM>& mf, IntVect const& nghost, Real time,
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

    FillPatchTwoLevels_doit(mf,nghost,time,cmf,ct,fmf,ft,
                            scomp,dcomp,ncomp,cgeom,fgeom,
                            cbc,cbccomp,fbc,fbccomp,ratio,mapper,bcs,bcscomp,
                            pre_interp,post_interp,index_space);
}

template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
EnableIf_t<IsFabArray<MF>::value>
FillPatchTwoLevels (Array<MF*, AMREX_SPACEDIM>& mf, Real time,
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

    FillPatchTwoLevels_doit(mf,mf.nGrowVect(),time,cmf,ct,fmf,ft,
                            scomp,dcomp,ncomp,cgeom,fgeom,
                            cbc,cbccomp,fbc,fbccomp,ratio,mapper,bcs,bcscomp,
                            pre_interp,post_interp,index_space);
}

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
                    const PostInterpHook& post_interp)
{
    FillPatchTwoLevels_doit(mf,nghost,time,cmf,ct,fmf,ft,
                            scomp,dcomp,ncomp,cgeom,fgeom,
                            cbc,cbccomp,fbc,fbccomp,ratio,mapper,bcs,bcscomp,
                            pre_interp,post_interp,&index_space);
}

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
                    const PostInterpHook& post_interp)
{
    FillPatchTwoLevels_doit(mf,mf.nGrowVect(),time,cmf,ct,fmf,ft,
                            scomp,dcomp,ncomp,cgeom,fgeom,
                            cbc,cbccomp,fbc,fbccomp,ratio,mapper,bcs,bcscomp,
                            pre_interp,post_interp,&index_space);
}
#endif

template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
EnableIf_t<IsFabArray<MF>::value>
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
EnableIf_t<IsFabArray<MF>::value>
InterpFromCoarseLevel (Array<MF*, AMREX_SPACEDIM>& mf, Real time,
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
EnableIf_t<IsFabArray<MF>::value>
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

    Box fdomain = fgeom.Domain();
    fdomain.convert(typ);

    Box fdomain_g(fdomain);
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

    mf_crse_patch.copy(cmf, scomp, 0, ncomp, cgeom.periodicity());

    cbc(mf_crse_patch, 0, ncomp, mf_crse_patch.nGrowVect(), time, cbccomp);

    int idummy1=0, idummy2=0;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        Vector<BCRec> bcr(ncomp);

        for (MFIter mfi(mf_crse_patch); mfi.isValid(); ++mfi)
        {
            FAB& sfab = mf_crse_patch[mfi];
            FAB& dfab = mf[mfi];
            Box dfab_bx = dfab.box();
            dfab_bx.grow(nghost-mf.nGrowVect());
            const Box& dbx = dfab_bx & fdomain_g;

            amrex::setBC(dbx,fdomain,bcscomp,0,ncomp,bcs,bcr);

            pre_interp(sfab, sfab.box(), 0, ncomp);

            mapper->interp(sfab,
                           0,
                           dfab,
                           dcomp,
                           ncomp,
                           dbx,
                           ratio,
                           cgeom,
                           fgeom,
                           bcr,
                           idummy1, idummy2, RunOn::Gpu);

            post_interp(dfab, dbx, dcomp, ncomp);
        }
    }

    fbc(mf, dcomp, ncomp, nghost, time, fbccomp);
}

template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
EnableIf_t<IsFabArray<MF>::value>
InterpFromCoarseLevel (Array<MF*, AMREX_SPACEDIM>& mf, IntVect const& nghost, Real time,
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

    const InterpolaterBoxCoarsener& coarsener = mapper->BoxCoarsener(ratio);
    const BoxArray& ba = mf[0]->boxArray();
    const DistributionMapping& dm = mf[0]->DistributionMap();

    // Change asserts.
    // If DM is different? Need to adjust or do all 3 seperately?

    // Confirm each boxArray is equal except for indexing?
    //     Note: if this is not true, have to setup all 3 directions separately.

    // const IndexType& typ = ba.ixType();
    // BL_ASSERT(typ == cmf.boxArray().ixType());

    // Create a cell-centered boxArray of the region to interp.
    // Convert this boxArray and domain as needed.
    Box fdomain = amrex::convert(fgeom.Domain(), IntVect::TheZeroVector());
    Box fdomain_g(fdomain);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (fgeom.isPeriodic(i)) {
            fdomain_g.grow(i,nghost[i]);
        }
    }

    BoxArray ba_crse_patch(ba.size());
    {  // TODO: later we might want to cache this
        for (int i = 0, N = ba.size(); i < N; ++i)
        {
            Box bx = amrex::convert(amrex::grow(ba[i], nghost), IntVect::TheZeroVector());
            bx &= fdomain_g;
            ba_crse_patch.set(i, coarsener.doit(bx));
        }
    }

    Array<MF, AMREX_SPACEDIM> mf_crse_patch;
    for (int i = 0; i<AMREX_SPACEDIM; ++i)
    {
        BoxArray ba_crse_idxed = amrex::convert(ba_crse_patch, mf[i]->boxArray().ixType());

#ifdef AMREX_USE_EB
        auto factory = makeEBFabFactory(cgeom, ba_crse_idxed, dm, {0,0,0}, EBSupport::basic);
        mf_crse_patch[i].define(ba_crse_idxed, dm, ncomp, 0, MFInfo(), *factory);
#else
        mf_crse_patch[i].define(ba_crse_idxed, dm, ncomp, 0);
#endif
        mf_set_domain_bndry(mf_crse_patch[i], cgeom);

        mf_crse_patch[i].copy(*(cmf[i]), scomp, 0, ncomp, cgeom.periodicity());

        cbc[i](mf_crse_patch[i], 0, ncomp, mf_crse_patch[i].nGrowVect(), time, cbccomp);
    }

    int idummy1=0, idummy2=0;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        Vector<Array<BCRec, AMREX_SPACEDIM> > bcr(ncomp);

        for (MFIter mfi(mf_crse_patch[0]); mfi.isValid(); ++mfi)
        {
            Array<FAB*, AMREX_SPACEDIM> sfab{ AMREX_D_DECL( &(mf_crse_patch[0][mfi]),
                                                            &(mf_crse_patch[1][mfi]),
                                                            &(mf_crse_patch[2][mfi]) )};
            Array<FAB*, AMREX_SPACEDIM> dfab{ AMREX_D_DECL( &((*mf[0])[mfi]),
                                                            &((*mf[1])[mfi]),
                                                            &((*mf[2])[mfi]) )};

            const Box& sbx_cc = amrex::convert(sfab[0]->box(), IntVect::TheZeroVector());
            Box dfab_cc = amrex::convert(dfab[0]->box(), IntVect::TheZeroVector());
            dfab_cc.grow(nghost-mf[0]->nGrowVect());
            const Box& dbx_cc = dfab_cc & fdomain_g;

            for (int d=0; d<AMREX_SPACEDIM; ++d)
            {
                Vector<BCRec> bcr_d(ncomp);

                amrex::setBC(dfab[d]->box(), amrex::convert(fdomain, dfab[d]->box().ixType()),
                             bcscomp,0,ncomp,bcs[d],bcr_d);

                for (int n=0; n<ncomp; ++n)
                    { bcr[n][d] = bcr_d[n]; }
            }

            pre_interp(sfab, sbx_cc, 0, ncomp);

            mapper->interp_arr(sfab, 0,  dfab, dcomp, ncomp, dbx_cc, ratio,
                               cgeom, fgeom, bcr, idummy1, idummy2, RunOn::Gpu);

            post_interp(dfab, dbx_cc, dcomp, ncomp);
        }
    }

    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
        fbc[d](*mf[d], dcomp, ncomp, nghost, time, fbccomp);
    }
}

}

#endif
