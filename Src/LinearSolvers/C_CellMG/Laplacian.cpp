
#include <winstd.H>
#include <Laplacian.H>
#include <LP_F.H>

Laplacian::Laplacian (const BndryData& bd,
                      Real             _h)
    :
    LinOp(bd,_h) {}

Laplacian::~Laplacian() {}

Real
Laplacian::norm (int nm, int level, const bool local)
{
  switch ( nm )
    {
    case 0:
      return 8.0/(h[level][0]*h[level][0]);
    }
  BoxLib::Error("Bad Laplacian::norm");
  return -1.0;
}

void
Laplacian::compFlux (D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
		     MultiFab& in, const BC_Mode& bc_mode)
{
    const int level    = 0;
    const int src_comp = 0;
    const int num_comp = 1;
    applyBC(in,src_comp,num_comp,level,bc_mode);
    const int nc = in.nComp();

    for (MFIter inmfi(in); inmfi.isValid(); ++inmfi)
    {
        const Box& vbx   = inmfi.validbox();
        FArrayBox& infab = in[inmfi];

        D_TERM(FArrayBox& xfab  = xflux[inmfi];,
               FArrayBox& yfab  = yflux[inmfi];,
               FArrayBox& zfab  = zflux[inmfi];);

        FORT_FLUX(infab.dataPtr(),
		  ARLIM(infab.loVect()), ARLIM(infab.hiVect()),
		  vbx.loVect(), vbx.hiVect(), &nc,
		  h[level],
		  xfab.dataPtr(),
		  ARLIM(xfab.loVect()), ARLIM(xfab.hiVect())
#if (BL_SPACEDIM >= 2)
		  ,yfab.dataPtr(),
		  ARLIM(yfab.loVect()), ARLIM(yfab.hiVect())
#endif
#if (BL_SPACEDIM == 3)
		  ,zfab.dataPtr(),
		  ARLIM(zfab.loVect()), ARLIM(zfab.hiVect())
#endif
		  );
    }
}

void
Laplacian::Fsmooth (MultiFab&       solnL,
                    const MultiFab& rhsL,
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

    const int nc = rhsL.nComp();

    for (MFIter solnLmfi(solnL); solnLmfi.isValid(); ++solnLmfi)
    {
        oitr.rewind();

        const int gn = solnLmfi.index();

        const LinOp::MaskTuple& mtuple = maskvals[level][gn];

        const Mask& m0 = *mtuple[oitr()]; oitr++;
        const Mask& m1 = *mtuple[oitr()]; oitr++;
        const Mask& m2 = *mtuple[oitr()]; oitr++;
        const Mask& m3 = *mtuple[oitr()]; oitr++;
#if (BL_SPACEDIM > 2 )
        const Mask& m4 = *mtuple[oitr()]; oitr++;
        const Mask& m5 = *mtuple[oitr()]; oitr++;
#endif
        const Box&       vbx     = solnLmfi.validbox();
        FArrayBox&       solnfab = solnL[gn];
        const FArrayBox& rhsfab  = rhsL[gn];
        const FArrayBox& f0fab   = f0[gn];
        const FArrayBox& f1fab   = f1[gn];
        const FArrayBox& f2fab   = f2[gn];
        const FArrayBox& f3fab   = f3[gn];
#if (BL_SPACEDIM == 3)
        const FArrayBox& f4fab   = f4[gn];
        const FArrayBox& f5fab   = f5[gn];
#endif

#if (BL_SPACEDIM == 2)
        FORT_GSRB(
            solnfab.dataPtr(), 
            ARLIM(solnfab.loVect()),ARLIM(solnfab.hiVect()),
            rhsfab.dataPtr(), 
            ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
            f0fab.dataPtr(), 
            ARLIM(f0fab.loVect()), ARLIM(f0fab.hiVect()),
            m0.dataPtr(), 
            ARLIM(m0.loVect()), ARLIM(m0.hiVect()),
            f1fab.dataPtr(), 
            ARLIM(f1fab.loVect()), ARLIM(f1fab.hiVect()),
            m1.dataPtr(), 
            ARLIM(m1.loVect()), ARLIM(m1.hiVect()),
            f2fab.dataPtr(), 
            ARLIM(f2fab.loVect()), ARLIM(f2fab.hiVect()),
            m2.dataPtr(), 
            ARLIM(m2.loVect()), ARLIM(m2.hiVect()),
            f3fab.dataPtr(), 
            ARLIM(f3fab.loVect()), ARLIM(f3fab.hiVect()),
            m3.dataPtr(), 
            ARLIM(m3.loVect()), ARLIM(m3.hiVect()),
            vbx.loVect(), vbx.hiVect(), &nc,
            h[level], &redBlackFlag);
#endif

#if (BL_SPACEDIM == 3)
        FORT_GSRB(
            solnfab.dataPtr(), 
            ARLIM(solnfab.loVect()),ARLIM(solnfab.hiVect()),
            rhsfab.dataPtr(), 
            ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
            f0fab.dataPtr(), 
            ARLIM(f0fab.loVect()), ARLIM(f0fab.hiVect()),
            m0.dataPtr(), 
            ARLIM(m0.loVect()), ARLIM(m0.hiVect()),
            f1fab.dataPtr(), 
            ARLIM(f1fab.loVect()), ARLIM(f1fab.hiVect()),
            m1.dataPtr(), 
            ARLIM(m1.loVect()), ARLIM(m1.hiVect()),
            f2fab.dataPtr(), 
            ARLIM(f2fab.loVect()), ARLIM(f2fab.hiVect()),
            m2.dataPtr(), 
            ARLIM(m2.loVect()), ARLIM(m2.hiVect()),
            f3fab.dataPtr(), 
            ARLIM(f3fab.loVect()), ARLIM(f3fab.hiVect()),
            m3.dataPtr(), 
            ARLIM(m3.loVect()), ARLIM(m3.hiVect()),
            f4fab.dataPtr(), 
            ARLIM(f4fab.loVect()), ARLIM(f4fab.hiVect()),
            m4.dataPtr(), 
            ARLIM(m4.loVect()), ARLIM(m4.hiVect()),
            f5fab.dataPtr(), 
            ARLIM(f5fab.loVect()), ARLIM(f5fab.hiVect()),
            m5.dataPtr(), 
            ARLIM(m5.loVect()), ARLIM(m5.hiVect()),
            vbx.loVect(), vbx.hiVect(), &nc,
            h[level], &redBlackFlag);
#endif
    }
}

void
Laplacian::Fsmooth_jacobi (MultiFab&       solnL,
                           const MultiFab& rhsL,
                           int            level)
{
}

void
Laplacian::Fapply (MultiFab&       y,
                   const MultiFab& x,
                   int             level)
{
    const int nc = y.nComp();

    for (MFIter ymfi(y); ymfi.isValid(); ++ymfi)
    {
        const Box&       vbx  = ymfi.validbox();
        FArrayBox&       yfab = y[ymfi];
        const FArrayBox& xfab = x[ymfi];

        FORT_ADOTX(yfab.dataPtr(), 
                   ARLIM(yfab.loVect()), ARLIM(yfab.hiVect()),
                   xfab.dataPtr(), 
                   ARLIM(xfab.loVect()), ARLIM(xfab.hiVect()),
                   vbx.loVect(), vbx.hiVect(), &nc,
                   h[level]);
    }
}
