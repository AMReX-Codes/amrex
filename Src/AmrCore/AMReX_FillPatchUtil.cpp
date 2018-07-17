#include <AMReX_Utility.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FillPatchUtil_F.H>
#include <cmath>
#include <limits>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex
{
    bool ProperlyNested (const IntVect& ratio, const IntVect& blocking_factor, int ngrow,
			 const IndexType& boxType, Interpolater* mapper)
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

    void FillPatchSingleLevel (MultiFab& mf, Real time, 
			       const Vector<MultiFab*>& smf, const Vector<Real>& stime,
			       int scomp, int dcomp, int ncomp,
			       const Geometry& geom, PhysBCFunctBase& physbcf)
    {
	BL_PROFILE("FillPatchSingleLevel");

	BL_ASSERT(scomp+ncomp <= smf[0]->nComp());
	BL_ASSERT(dcomp+ncomp <= mf.nComp());
	BL_ASSERT(smf.size() == stime.size());
	BL_ASSERT(smf.size() != 0);

	if (smf.size() == 1) 
	{
	    mf.copy(*smf[0], scomp, dcomp, ncomp, 0, mf.nGrow(), geom.periodicity());
	} 
	else if (smf.size() == 2) 
	{
	    BL_ASSERT(smf[0]->boxArray() == smf[1]->boxArray());
	    MultiFab raii;
	    MultiFab * dmf;
	    int destcomp;
	    bool sameba;
	    if (mf.boxArray() == smf[0]->boxArray()) {
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

#ifdef _OPENMP
#pragma omp parallel 
#endif
	    for (MFIter mfi(*dmf,true); mfi.isValid(); ++mfi)
	    {
		const Box& bx = mfi.tilebox();
		(*dmf)[mfi].linInterp((*smf[0])[mfi],
				      scomp,
				      (*smf[1])[mfi],
				      scomp,
				      stime[0],
				      stime[1],
				      time,
				      bx,
				      destcomp,
				      ncomp);
	    }
	    
	    if (sameba)
	    {
		// Note that when sameba is true mf's BoxArray is nonoverlapping.
		// So FillBoundary is safe.
		mf.FillBoundary(dcomp,ncomp,geom.periodicity());
	    }
	    else
	    {
		int src_ngrow = 0;
		int dst_ngrow = mf.nGrow();

		mf.copy(*dmf, 0, dcomp, ncomp, src_ngrow, dst_ngrow, geom.periodicity());
	    }
	}
	else {
	    amrex::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
	}

	physbcf.FillBoundary(mf, dcomp, ncomp, time);
    }

    void FillPatchTwoLevels (MultiFab& mf, Real time,
			     const Vector<MultiFab*>& cmf, const Vector<Real>& ct,
			     const Vector<MultiFab*>& fmf, const Vector<Real>& ft,
			     int scomp, int dcomp, int ncomp,
			     const Geometry& cgeom, const Geometry& fgeom, 
			     PhysBCFunctBase& cbc, PhysBCFunctBase& fbc,
			     const IntVect& ratio, 
			     Interpolater* mapper, const BCRec& bcs,
                             const InterpHook& pre_interp,
                             const InterpHook& post_interp)
    {
        Vector<BCRec> bcs_array(1,BCRec(bcs.lo(),bcs.hi()));

        FillPatchTwoLevels(mf,time,cmf,ct,fmf,ft,scomp,dcomp,ncomp,cgeom,fgeom,
                           cbc,fbc,ratio,mapper,bcs_array,pre_interp,post_interp);
    }


    void FillPatchTwoLevels (MultiFab& mf, Real time,
			     const Vector<MultiFab*>& cmf, const Vector<Real>& ct,
			     const Vector<MultiFab*>& fmf, const Vector<Real>& ft,
			     int scomp, int dcomp, int ncomp,
			     const Geometry& cgeom, const Geometry& fgeom, 
			     PhysBCFunctBase& cbc, PhysBCFunctBase& fbc,
			     const IntVect& ratio, 
			     Interpolater* mapper, const Vector<BCRec>& bcs,
                             const InterpHook& pre_interp,
                             const InterpHook& post_interp)
    {
	BL_PROFILE("FillPatchTwoLevels");

	int ngrow = mf.nGrow();
	    
	if (ngrow > 0 || mf.getBDKey() != fmf[0]->getBDKey()) 
	{
	    const InterpolaterBoxCoarsener& coarsener = mapper->BoxCoarsener(ratio);
	    
	    Box fdomain = fgeom.Domain();
	    fdomain.convert(mf.boxArray().ixType());
	    Box fdomain_g(fdomain);
	    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
		if (fgeom.isPeriodic(i)) {
		    fdomain_g.grow(i,ngrow);
		}
	    }
	    
	    const FabArrayBase::FPinfo& fpc = FabArrayBase::TheFPinfo(*fmf[0], mf, fdomain_g,
                                                                      IntVect(ngrow),
                                                                      coarsener, 
                                                                      amrex::coarsen(fgeom.Domain(),ratio));

	    if ( ! fpc.ba_crse_patch.empty())
	    {
		MultiFab mf_crse_patch(fpc.ba_crse_patch, fpc.dm_crse_patch, ncomp, 0, MFInfo(),
                                       *fpc.fact_crse_patch);
		
                mf_crse_patch.setDomainBndry(std::numeric_limits<Real>::quiet_NaN(), cgeom);

		FillPatchSingleLevel(mf_crse_patch, time, cmf, ct, scomp, 0, ncomp, cgeom, cbc);

		int idummy1=0, idummy2=0;
		bool cc = fpc.ba_crse_patch.ixType().cellCentered();
                ignore_unused(cc);
#ifdef _OPENMP
#pragma omp parallel if (cc)
#endif
		for (MFIter mfi(mf_crse_patch); mfi.isValid(); ++mfi)
		{
                    FArrayBox& sfab = mf_crse_patch[mfi];
		    int li = mfi.LocalIndex();
		    int gi = fpc.dst_idxs[li];	
                    FArrayBox& dfab = mf[gi];
		    const Box& dbx = fpc.dst_boxes[li];
		    
		    Vector<BCRec> bcr(ncomp);
		    amrex::setBC(dbx,fdomain,scomp,0,ncomp,bcs,bcr);

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
				   idummy1, idummy2);

                    post_interp(dfab, dbx, dcomp, ncomp);
		}
	    }
	}

	FillPatchSingleLevel(mf, time, fmf, ft, scomp, dcomp, ncomp, fgeom, fbc);
    }

    void InterpFromCoarseLevel (MultiFab& mf, Real time, const MultiFab& cmf, 
				int scomp, int dcomp, int ncomp,
				const Geometry& cgeom, const Geometry& fgeom, 
				PhysBCFunctBase& cbc, PhysBCFunctBase& fbc, const IntVect& ratio, 
				Interpolater* mapper, const BCRec& bcs,
                                const InterpHook& pre_interp,
                                const InterpHook& post_interp)
    {

        Vector<BCRec> bcs_array(1,BCRec(bcs.lo(),bcs.hi()));
        InterpFromCoarseLevel(mf,time,cmf,scomp,dcomp,ncomp,cgeom,fgeom,
                              cbc,fbc,ratio,mapper,bcs_array,
                              pre_interp, post_interp);

    }


    void InterpFromCoarseLevel (MultiFab& mf, Real time, const MultiFab& cmf, 
				int scomp, int dcomp, int ncomp,
				const Geometry& cgeom, const Geometry& fgeom, 
				PhysBCFunctBase& cbc, PhysBCFunctBase& fbc, const IntVect& ratio, 
				Interpolater* mapper, const Vector<BCRec>& bcs,
                                const InterpHook& pre_interp,
                                const InterpHook& post_interp)
    {
	const InterpolaterBoxCoarsener& coarsener = mapper->BoxCoarsener(ratio);

	const BoxArray& ba = mf.boxArray();
	const DistributionMapping& dm = mf.DistributionMap();
	int ngrow = mf.nGrow();

	const IndexType& typ = ba.ixType();

	BL_ASSERT(typ == cmf.boxArray().ixType());

	Box fdomain = fgeom.Domain();
	fdomain.convert(typ);

	Box fdomain_g(fdomain);
	for (int i = 0; i < AMREX_SPACEDIM; ++i) {
	    if (fgeom.isPeriodic(i)) {
		fdomain_g.grow(i,ngrow);
	    }
	}

	BoxArray ba_crse_patch(ba.size());
	{  // TODO: later we might want to cache this
	    for (int i = 0, N = ba.size(); i < N; ++i)
	    {
		Box bx = amrex::convert(amrex::grow(ba[i],ngrow), typ);
		bx &= fdomain_g;
		ba_crse_patch.set(i, coarsener.doit(bx));
	    }
	}

#ifdef AMREX_USE_EB
        auto factory = makeEBFabFactory(cgeom, ba_crse_patch, dm, {0,0,0}, EBSupport::basic);
	MultiFab mf_crse_patch(ba_crse_patch, dm, ncomp, 0, MFInfo(), *factory);
#else
	MultiFab mf_crse_patch(ba_crse_patch, dm, ncomp, 0);
#endif

        mf_crse_patch.setDomainBndry(std::numeric_limits<Real>::quiet_NaN(), cgeom);

	mf_crse_patch.copy(cmf, scomp, 0, ncomp, cgeom.periodicity());

	cbc.FillBoundary(mf_crse_patch, 0, ncomp, time);

	int idummy1=0, idummy2=0;

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
	    Vector<BCRec> bcr(ncomp);

            for (MFIter mfi(mf_crse_patch); mfi.isValid(); ++mfi)
            {
                FArrayBox& sfab = mf_crse_patch[mfi];
                FArrayBox& dfab = mf[mfi];
                const Box& dbx = dfab.box() & fdomain_g;
                
                amrex::setBC(dbx,fdomain,scomp,0,ncomp,bcs,bcr);

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
                               idummy1, idummy2);

                post_interp(dfab, dbx, dcomp, ncomp);
            }
	}

	fbc.FillBoundary(mf, dcomp, ncomp, time);
    }

    // B fields are assumed to be on staggered grids.
    void InterpCrseFineBndryEMfield (InterpEM_t interp_type,
                                     const std::array<MultiFab,AMREX_SPACEDIM>& crse,
                                     std::array<MultiFab,AMREX_SPACEDIM>& fine,
                                     const Geometry& cgeom, const Geometry& fgeom,
                                     int ref_ratio)
    {
        BL_ASSERT(ref_ratio == 2);

        int ngrow = fine[0].nGrow();
        for (int idim = 1; idim < AMREX_SPACEDIM; ++idim) {
            BL_ASSERT(ngrow == fine[idim].nGrow());
        }

        if (ngrow == 0) return;

        bool include_periodic = true;
        bool include_physbndry = false;
        const auto& cfinfo = FabArrayBase::TheCFinfo(fine[0], fgeom, IntVect(ngrow),
                                                     include_periodic, include_physbndry);

        if (! cfinfo.ba_cfb.empty())
        {
            std::array<MultiFab, AMREX_SPACEDIM> cmf;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const BoxArray& fine_ba = fine[idim].boxArray();
                BoxArray fba = cfinfo.ba_cfb;
                fba.convert(fine_ba.ixType());
                BoxArray cba = fba;
                cba.coarsen(ref_ratio);
                const DistributionMapping& dm = cfinfo.dm_cfb;

#ifdef AMREX_USE_EB
                amrex::Abort("InterpCrseFineBndryEMfield: EB is allowed");
#endif

                cmf[idim].define(cba, dm, 1, 1, MFInfo(), crse[0].Factory());

                cmf[idim].copy(crse[idim], 0, 0, 1, 0, 1, cgeom.periodicity());
            }

            const Real* dx = cgeom.CellSize();

            const int use_limiter = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
            {
                std::array<FArrayBox,AMREX_SPACEDIM> bfab;
                for (MFIter mfi(cmf[0]); mfi.isValid(); ++mfi)
                {
                    const int fi = cfinfo.fine_grid_idx[mfi.LocalIndex()];

                    Box ccbx = amrex::grow(fine[0].boxArray()[fi], ngrow);
                    ccbx.enclosedCells();
                    ccbx.coarsen(ref_ratio).refine(ref_ratio);  // so that ccbx is coarsenable

                    const FArrayBox& cxfab = cmf[0][mfi];
                    bfab[0].resize(amrex::convert(ccbx,fine[0].ixType()));
#if (AMREX_SPACEDIM > 1)
                    const FArrayBox& cyfab = cmf[1][mfi];
                    bfab[1].resize(amrex::convert(ccbx,fine[1].ixType()));
#endif
#if (AMREX_SPACEDIM > 2)
                    const FArrayBox& czfab = cmf[2][mfi];
                    bfab[2].resize(amrex::convert(ccbx,fine[2].ixType()));
#endif
                    // interpolate from cmf to fmf
                    if (interp_type == InterpB)
                    {
                        amrex_interp_div_free_bfield(BL_TO_FORTRAN_BOX(ccbx),
                                                     AMREX_D_DECL(BL_TO_FORTRAN_ANYD(bfab[0]),
                                                                  BL_TO_FORTRAN_ANYD(bfab[1]),
                                                                  BL_TO_FORTRAN_ANYD(bfab[2])),
                                                     AMREX_D_DECL(BL_TO_FORTRAN_ANYD(cxfab),
                                                                  BL_TO_FORTRAN_ANYD(cyfab),
                                                                  BL_TO_FORTRAN_ANYD(czfab)),
                                                     dx, &ref_ratio, &use_limiter);
                    }
                    else if (interp_type == InterpE)
                    {
                        amrex_interp_efield(BL_TO_FORTRAN_BOX(ccbx),
                                            AMREX_D_DECL(BL_TO_FORTRAN_ANYD(bfab[0]),
                                                         BL_TO_FORTRAN_ANYD(bfab[1]),
                                                         BL_TO_FORTRAN_ANYD(bfab[2])),
                                            AMREX_D_DECL(BL_TO_FORTRAN_ANYD(cxfab),
                                                         BL_TO_FORTRAN_ANYD(cyfab),
                                                         BL_TO_FORTRAN_ANYD(czfab)),
                                            &ref_ratio, &use_limiter);
                    }
                    else
                    {
                        amrex::Abort("InterpCrseFineBndryEMfield: unknown interp_type");
                    }
                    
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                    {
                        const BoxArray& fine_ba = fine[idim].boxArray();
                        const Box& fine_valid_box = fine_ba[fi];
                        Box b = bfab[idim].box();
                        b &= fine_valid_box;
                        const BoxList& diff = amrex::boxDiff(b, fine_valid_box); // skip valid cells
                        FArrayBox& fine_fab = fine[idim][fi];
                        for (const auto& x : diff)
                        {
                            fine_fab.copy(bfab[idim], x, 0, x, 0, 1);
                        }
                    }
                }
            }
        }
    }
}
