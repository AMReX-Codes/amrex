#include <algorithm>
#include <AMReX_ABec2.H>
#include <AMReX_ABec2_F.H>
#include <AMReX_ParallelDescriptor.H>

#include <AMReX_LO_BCTYPES.H>
#include <AMReX_LO_F.H>

#include <AMReX_BLFort.H>

namespace amrex {

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

        BL_ASSERT(level<undrrelxr.size());

        const BndryData::RealTuple&      bdl = bgb->bndryLocs(gn);
        const Vector< Vector<BoundCond> >& bdc = bgb->bndryConds(gn);

        for (OrientationIter oitr; oitr; ++oitr)
        {
            const Orientation o = oitr();

            FabSet&       f   = undrrelxr[level][o];
            int           cdr = o;
            const Mask&   m   = local ? lmaskvals[level][o][mfi] : maskvals[level][o][mfi];
            Real          bcl = bdl[o];
            BL_ASSERT(bdc[o].size()>bndry_comp);
            int           bct = bdc[o][bndry_comp];

            const Box&       vbx   = ba[gn];

            BL_ASSERT(f.size()>gn);
            FArrayBox&       ffab  = f[gn];

            amrex_ab2_bndryrlx
              ( vbx.loVect(), vbx.hiVect(),
                BL_TO_FORTRAN(ffab),
                BL_TO_FORTRAN(m),
                &cdr, &bct, &bcl, &maxorder, h[level].data());
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

  const FabSet& f0 = undrrelxr[level][oitr()]; oitr++;
  const FabSet& f1 = undrrelxr[level][oitr()]; oitr++;
  const FabSet& f2 = undrrelxr[level][oitr()]; oitr++;
  const FabSet& f3 = undrrelxr[level][oitr()]; oitr++;
#if (BL_SPACEDIM > 2)
  const FabSet& f4 = undrrelxr[level][oitr()]; oitr++;
  const FabSet& f5 = undrrelxr[level][oitr()]; oitr++;
#endif    
  const MultiFab& a = aCoefficients(level);

  AMREX_D_TERM(const MultiFab& bX = bCoefficients(0,level);,
         const MultiFab& bY = bCoefficients(1,level);,
         const MultiFab& bZ = bCoefficients(2,level););

  oitr.rewind();
  const MultiMask& mm0 = maskvals[level][oitr()]; oitr++;
  const MultiMask& mm1 = maskvals[level][oitr()]; oitr++;
  const MultiMask& mm2 = maskvals[level][oitr()]; oitr++;
  const MultiMask& mm3 = maskvals[level][oitr()]; oitr++;
#if (BL_SPACEDIM > 2)
  const MultiMask& mm4 = maskvals[level][oitr()]; oitr++;
  const MultiMask& mm5 = maskvals[level][oitr()]; oitr++;
#endif

  const int nc = 1;

  Real alpha = get_alpha();
  Real beta = get_beta();

  const bool tiling = true;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter solnLmfi(solnL,tiling); solnLmfi.isValid(); ++solnLmfi)
  {
    const Mask& m0 = mm0[solnLmfi];
    const Mask& m1 = mm1[solnLmfi];
    const Mask& m2 = mm2[solnLmfi];
    const Mask& m3 = mm3[solnLmfi];
#if (BL_SPACEDIM > 2)
    const Mask& m4 = mm4[solnLmfi];
    const Mask& m5 = mm5[solnLmfi];
#endif

    const Box&       tbx     = solnLmfi.tilebox();
    const Box&       vbx     = solnLmfi.validbox();
    FArrayBox&       solnfab = solnL[solnLmfi];
    const FArrayBox& resfab  = resL[solnLmfi];
    const FArrayBox& afab    = a[solnLmfi];

    AMREX_D_TERM(const FArrayBox& bxfab = bX[solnLmfi];,
           const FArrayBox& byfab = bY[solnLmfi];,
           const FArrayBox& bzfab = bZ[solnLmfi];);

    const FArrayBox& f0fab = f0[solnLmfi];
    const FArrayBox& f1fab = f1[solnLmfi];
    const FArrayBox& f2fab = f2[solnLmfi];
    const FArrayBox& f3fab = f3[solnLmfi];
#if (BL_SPACEDIM == 3)
    const FArrayBox& f4fab = f4[solnLmfi];
    const FArrayBox& f5fab = f5[solnLmfi];
#endif

    amrex_ab2_gsrb
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
        &nc, h[level].data(), &redBlackFlag);

  }
}

void
ABec2::Fsmooth_jacobi (MultiFab&       solnL,
		       const MultiFab& resL,
		       int             level)
{
  OrientationIter oitr;

  const FabSet& f0 = undrrelxr[level][oitr()]; oitr++;
  const FabSet& f1 = undrrelxr[level][oitr()]; oitr++;
  const FabSet& f2 = undrrelxr[level][oitr()]; oitr++;
  const FabSet& f3 = undrrelxr[level][oitr()]; oitr++;
#if (BL_SPACEDIM > 2)
  const FabSet& f4 = undrrelxr[level][oitr()]; oitr++;
  const FabSet& f5 = undrrelxr[level][oitr()]; oitr++;
#endif    
  const MultiFab& a = aCoefficients(level);

  AMREX_D_TERM(const MultiFab& bX = bCoefficients(0,level);,
         const MultiFab& bY = bCoefficients(1,level);,
         const MultiFab& bZ = bCoefficients(2,level););

  oitr.rewind();
  const MultiMask& mm0 = maskvals[level][oitr()]; oitr++;
  const MultiMask& mm1 = maskvals[level][oitr()]; oitr++;
  const MultiMask& mm2 = maskvals[level][oitr()]; oitr++;
  const MultiMask& mm3 = maskvals[level][oitr()]; oitr++;
#if (BL_SPACEDIM > 2)
  const MultiMask& mm4 = maskvals[level][oitr()]; oitr++;
  const MultiMask& mm5 = maskvals[level][oitr()]; oitr++;
#endif

  //const int nc = solnL.nComp(); // FIXME: This LinOp only really supports single-component
  const int nc = 1;
  Real alpha = get_alpha();
  Real beta = get_beta();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter solnLmfi(solnL); solnLmfi.isValid(); ++solnLmfi)
  {
    const Mask& m0 = mm0[solnLmfi];
    const Mask& m1 = mm1[solnLmfi];
    const Mask& m2 = mm2[solnLmfi];
    const Mask& m3 = mm3[solnLmfi];
#if (BL_SPACEDIM > 2)
    const Mask& m4 = mm4[solnLmfi];
    const Mask& m5 = mm5[solnLmfi];
#endif

    const Box&       vbx     = solnLmfi.validbox();
    FArrayBox&       solnfab = solnL[solnLmfi];
    const FArrayBox& resfab  = resL[solnLmfi];
    const FArrayBox& afab    = a[solnLmfi];

    AMREX_D_TERM(const FArrayBox& bxfab = bX[solnLmfi];,
           const FArrayBox& byfab = bY[solnLmfi];,
           const FArrayBox& bzfab = bZ[solnLmfi];);

    const FArrayBox& f0fab = f0[solnLmfi];
    const FArrayBox& f1fab = f1[solnLmfi];
    const FArrayBox& f2fab = f2[solnLmfi];
    const FArrayBox& f3fab = f3[solnLmfi];
#if (BL_SPACEDIM == 3)
    const FArrayBox& f4fab = f4[solnLmfi];
    const FArrayBox& f5fab = f5[solnLmfi];
#endif

    amrex_ab2_jacobi
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
        &nc, h[level].data());

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
    amrex::Abort("Shouldnt be here");
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
    amrex::Abort("Shouldnt be here");
  }
}

}

