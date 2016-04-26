#include <winstd.H>
#include <algorithm>
#include <ABec2.H>
#include <ABec2_F.H>
#include <ParallelDescriptor.H>

#include <LO_BCTYPES.H>
#include <LO_F.H>

#include <BLFort.H>

void
ABec2::altApplyBC (int  level,
                   bool local)
{
    //
    // No coarsened boundary values, cannot apply inhomog at lev>0.
    //
    BL_ASSERT(level < numLevels());

    prepareForLevel(level);
    //
    // Fill boundary cells.
    //
    const MultiFab& a = aCoefficients(level);

    int bndry_comp = 0;
    const BoxArray& ba = boxArray(level);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(a); mfi.isValid(); ++mfi)
    {
        const int gn = mfi.index();

        BL_ASSERT(level<maskvals.size() && maskvals[level].find(gn)!=maskvals[level].end());
        BL_ASSERT(level<lmaskvals.size() && lmaskvals[level].find(gn)!=lmaskvals[level].end());
        BL_ASSERT(level<undrrelxr.size());

        const MaskTuple&                 ma  =  maskvals[level][gn];
        const MaskTuple&                 lma = lmaskvals[level][gn];
        const BndryData::RealTuple&      bdl = bgb->bndryLocs(gn);
        const Array< Array<BoundCond> >& bdc = bgb->bndryConds(gn);

        for (OrientationIter oitr; oitr; ++oitr)
        {
            const Orientation o = oitr();

            FabSet&       f   = (*undrrelxr[level])[o];
            int           cdr = o;
            const Mask&   m   = local ? (*lma[o]) : (*ma[o]);
            Real          bcl = bdl[o];
            BL_ASSERT(bdc[o].size()>bndry_comp);
            int           bct = bdc[o][bndry_comp];

            const Box&       vbx   = ba[gn];

            BL_ASSERT(f.size()>gn);
            FArrayBox&       ffab  = f[gn];

            BL_FORT_PROC_CALL(AB2_BNDRYRLX, ab2_bndryrlx)
              ( vbx.loVect(), vbx.hiVect(),
                BL_TO_FORTRAN(ffab),
                BL_TO_FORTRAN(m),
                &cdr, &bct, &bcl, &maxorder, h[level]);
        }
    }
}



void
ABec2::Fsmooth (MultiFab&       solnL,
		const MultiFab& resL,
		int             level,
		int             redBlackFlag)
{
  OrientationIter oitr;

  const FabSet& f0 = (*undrrelxr[level])[oitr()]; oitr++;
  const FabSet& f1 = (*undrrelxr[level])[oitr()]; oitr++;
  const FabSet& f2 = (*undrrelxr[level])[oitr()]; oitr++;
  const FabSet& f3 = (*undrrelxr[level])[oitr()]; oitr++;
#if (BL_SPACEDIM > 2)
  const FabSet& f4 = (*undrrelxr[level])[oitr()]; oitr++;
  const FabSet& f5 = (*undrrelxr[level])[oitr()]; oitr++;
#endif    
  const MultiFab& a = aCoefficients(level);

  D_TERM(const MultiFab& bX = bCoefficients(0,level);,
         const MultiFab& bY = bCoefficients(1,level);,
         const MultiFab& bZ = bCoefficients(2,level););

  const int nc = 1;

  Real alpha = get_alpha();
  Real beta = get_beta();

  const bool tiling = true;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter solnLmfi(solnL,tiling); solnLmfi.isValid(); ++solnLmfi)
  {
    OrientationIter oitr;

    const int gn = solnLmfi.index();

    const LinOp::MaskTuple& mtuple = maskvals[level][gn];

    const Mask& m0 = *mtuple[oitr()]; oitr++;
    const Mask& m1 = *mtuple[oitr()]; oitr++;
    const Mask& m2 = *mtuple[oitr()]; oitr++;
    const Mask& m3 = *mtuple[oitr()]; oitr++;
#if (BL_SPACEDIM == 3)
    const Mask& m4 = *mtuple[oitr()]; oitr++;
    const Mask& m5 = *mtuple[oitr()]; oitr++;
#endif
    const Box&       tbx     = solnLmfi.tilebox();
    const Box&       vbx     = solnLmfi.validbox();
    FArrayBox&       solnfab = solnL[gn];
    const FArrayBox& resfab  = resL[gn];
    const FArrayBox& afab    = a[gn];

    D_TERM(const FArrayBox& bxfab = bX[gn];,
           const FArrayBox& byfab = bY[gn];,
           const FArrayBox& bzfab = bZ[gn];);

    const FArrayBox& f0fab = f0[gn];
    const FArrayBox& f1fab = f1[gn];
    const FArrayBox& f2fab = f2[gn];
    const FArrayBox& f3fab = f3[gn];
#if (BL_SPACEDIM == 3)
    const FArrayBox& f4fab = f4[gn];
    const FArrayBox& f5fab = f5[gn];
#endif

    BL_FORT_PROC_CALL(AB2_GSRB, ab2_gsrb)
      ( tbx.loVect(), tbx.hiVect(),
        vbx.loVect(), vbx.hiVect(),
        BL_TO_FORTRAN(solnfab),
        BL_TO_FORTRAN(resfab),
        &alpha, &beta,
        BL_TO_FORTRAN(afab),
        BL_TO_FORTRAN(bxfab),
	BL_TO_FORTRAN(byfab),
#if (BL_SPACEDIM == 3)
	BL_TO_FORTRAN(bzfab),
#endif
        BL_TO_FORTRAN(f0fab), BL_TO_FORTRAN(m0),
        BL_TO_FORTRAN(f1fab), BL_TO_FORTRAN(m1),
        BL_TO_FORTRAN(f2fab), BL_TO_FORTRAN(m2),
        BL_TO_FORTRAN(f3fab), BL_TO_FORTRAN(m3),
#if (BL_SPACEDIM == 3)
        BL_TO_FORTRAN(f4fab), BL_TO_FORTRAN(m4),
        BL_TO_FORTRAN(f5fab), BL_TO_FORTRAN(m5),
#endif
        &nc, h[level], &redBlackFlag);

  }
}

void
ABec2::Fsmooth_jacobi (MultiFab&       solnL,
		       const MultiFab& resL,
		       int             level)
{
  OrientationIter oitr;

  const FabSet& f0 = (*undrrelxr[level])[oitr()]; oitr++;
  const FabSet& f1 = (*undrrelxr[level])[oitr()]; oitr++;
  const FabSet& f2 = (*undrrelxr[level])[oitr()]; oitr++;
  const FabSet& f3 = (*undrrelxr[level])[oitr()]; oitr++;
#if (BL_SPACEDIM > 2)
  const FabSet& f4 = (*undrrelxr[level])[oitr()]; oitr++;
  const FabSet& f5 = (*undrrelxr[level])[oitr()]; oitr++;
#endif    
  const MultiFab& a = aCoefficients(level);

  D_TERM(const MultiFab& bX = bCoefficients(0,level);,
         const MultiFab& bY = bCoefficients(1,level);,
         const MultiFab& bZ = bCoefficients(2,level););

  //const int nc = solnL.nComp(); // FIXME: This LinOp only really supports single-component
  const int nc = 1;
  Real alpha = get_alpha();
  Real beta = get_beta();

  for (MFIter solnLmfi(solnL); solnLmfi.isValid(); ++solnLmfi)
  {
    oitr.rewind();

    const int gn = solnLmfi.index();

    const LinOp::MaskTuple& mtuple = maskvals[level][gn];

    const Mask& m0 = *mtuple[oitr()]; oitr++;
    const Mask& m1 = *mtuple[oitr()]; oitr++;
    const Mask& m2 = *mtuple[oitr()]; oitr++;
    const Mask& m3 = *mtuple[oitr()]; oitr++;
#if (BL_SPACEDIM == 3)
    const Mask& m4 = *mtuple[oitr()]; oitr++;
    const Mask& m5 = *mtuple[oitr()]; oitr++;
#endif
    const Box&       vbx     = solnLmfi.validbox();
    FArrayBox&       solnfab = solnL[gn];
    const FArrayBox& resfab  = resL[gn];
    const FArrayBox& afab    = a[gn];

    D_TERM(const FArrayBox& bxfab = bX[gn];,
           const FArrayBox& byfab = bY[gn];,
           const FArrayBox& bzfab = bZ[gn];);

    const FArrayBox& f0fab = f0[gn];
    const FArrayBox& f1fab = f1[gn];
    const FArrayBox& f2fab = f2[gn];
    const FArrayBox& f3fab = f3[gn];
#if (BL_SPACEDIM == 3)
    const FArrayBox& f4fab = f4[gn];
    const FArrayBox& f5fab = f5[gn];
#endif

    BL_FORT_PROC_CALL(AB2_JACOBI, ab2_jacobi)
      ( vbx.loVect(), vbx.hiVect(),
        BL_TO_FORTRAN(solnfab),
        BL_TO_FORTRAN(resfab),
        &alpha, &beta,
        BL_TO_FORTRAN(afab),
        BL_TO_FORTRAN(bxfab),
	BL_TO_FORTRAN(byfab),
#if (BL_SPACEDIM == 3)
	BL_TO_FORTRAN(bzfab),
#endif
        BL_TO_FORTRAN(f0fab), BL_TO_FORTRAN(m0),
        BL_TO_FORTRAN(f1fab), BL_TO_FORTRAN(m1),
        BL_TO_FORTRAN(f2fab), BL_TO_FORTRAN(m2),
        BL_TO_FORTRAN(f3fab), BL_TO_FORTRAN(m3),
#if (BL_SPACEDIM == 3)
        BL_TO_FORTRAN(f4fab), BL_TO_FORTRAN(m4),
        BL_TO_FORTRAN(f5fab), BL_TO_FORTRAN(m5),
#endif
        &nc, h[level]);

  }
}

void
ABec2::smooth (MultiFab&       solnL,
               const MultiFab& rhsL,
               int             level,
               LinOp::BC_Mode  bc_mode)
{
  if (level > 0)
  {
    for (int redBlackFlag = 0; redBlackFlag < 2; redBlackFlag++)
    {
      applyBC(solnL, 0, 1, level, bc_mode);
      ABecLaplacian::Fsmooth(solnL, rhsL, level, redBlackFlag);
    }
  }
  else
  {
    BoxLib::Abort("Shouldnt be here");
  }
}

void
ABec2::altSmooth (MultiFab&       solnL,
                  const MultiFab& resL,
                  int             level,
                  int             redBlackFlag)
{
  bool local = false;
  altApplyBC(level,local);
  Fsmooth(solnL, resL, level, redBlackFlag);
}

void
ABec2::jacobi_smooth (MultiFab&       solnL,
                      const MultiFab& rhsL,
                      int             level,
                      LinOp::BC_Mode  bc_mode)
{
  if (level > 0) {
    for (int redBlackFlag = 0; redBlackFlag < 2; redBlackFlag++)
    {
      applyBC(solnL, 0, 1, level, bc_mode);
      ABecLaplacian::Fsmooth_jacobi(solnL, rhsL, level);
    }
  }
  else
  {
    BoxLib::Abort("Shouldnt be here");
  }
}

