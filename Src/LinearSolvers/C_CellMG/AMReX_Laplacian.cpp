
#include <AMReX_Laplacian.H>
#include <AMReX_LP_F.H>

namespace amrex {

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
  amrex::Error("Bad Laplacian::norm");
  return -1.0;
}

void
Laplacian::compFlux (AMREX_D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
		     MultiFab& in, const BC_Mode& bc_mode,
		     int src_comp, int dst_comp, int num_comp, int bnd_comp)
{
    BL_PROFILE("Laplacian::compFlux()");

    const int level    = 0;
    applyBC(in,src_comp,num_comp,level,bc_mode,bnd_comp);

    const bool tiling = true;
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter inmfi(in,tiling); inmfi.isValid(); ++inmfi)
    {
        AMREX_D_TERM(const Box& xbx   = inmfi.nodaltilebox(0);,
	       const Box& ybx   = inmfi.nodaltilebox(1);,
	       const Box& zbx   = inmfi.nodaltilebox(2););

        FArrayBox& infab = in[inmfi];

        AMREX_D_TERM(FArrayBox& xfab  = xflux[inmfi];,
               FArrayBox& yfab  = yflux[inmfi];,
               FArrayBox& zfab  = zflux[inmfi];);

        amrex_lp_flux(infab.dataPtr(src_comp),
		  ARLIM(infab.loVect()), ARLIM(infab.hiVect()),
		  xbx.loVect(), xbx.hiVect(), 
#if (BL_SPACEDIM >= 2)
		  ybx.loVect(), ybx.hiVect(), 
#if (BL_SPACEDIM == 3)
		  zbx.loVect(), zbx.hiVect(), 
#endif
#endif
	          &num_comp,
		  h[level].data(),
		  xfab.dataPtr(dst_comp),
		  ARLIM(xfab.loVect()), ARLIM(xfab.hiVect())
#if (BL_SPACEDIM >= 2)
		  ,yfab.dataPtr(dst_comp),
		  ARLIM(yfab.loVect()), ARLIM(yfab.hiVect())
#endif
#if (BL_SPACEDIM == 3)
		  ,zfab.dataPtr(dst_comp),
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
    BL_PROFILE("Laplacian::Fsmooth()");

    OrientationIter oitr;

    const FabSet& f0 = undrrelxr[level][oitr()]; oitr++;
    const FabSet& f1 = undrrelxr[level][oitr()]; oitr++;
    const FabSet& f2 = undrrelxr[level][oitr()]; oitr++;
    const FabSet& f3 = undrrelxr[level][oitr()]; oitr++;
#if (BL_SPACEDIM > 2)
    const FabSet& f4 = undrrelxr[level][oitr()]; oitr++;
    const FabSet& f5 = undrrelxr[level][oitr()]; oitr++;
#endif

    oitr.rewind();
    const MultiMask& mm0 = maskvals[level][oitr()]; oitr++;
    const MultiMask& mm1 = maskvals[level][oitr()]; oitr++;
    const MultiMask& mm2 = maskvals[level][oitr()]; oitr++;
    const MultiMask& mm3 = maskvals[level][oitr()]; oitr++;
#if (BL_SPACEDIM > 2)
    const MultiMask& mm4 = maskvals[level][oitr()]; oitr++;
    const MultiMask& mm5 = maskvals[level][oitr()]; oitr++;
#endif

    const int nc = rhsL.nComp();

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
        const FArrayBox& rhsfab  = rhsL[solnLmfi];
        const FArrayBox& f0fab   = f0[solnLmfi];
        const FArrayBox& f1fab   = f1[solnLmfi];
        const FArrayBox& f2fab   = f2[solnLmfi];
        const FArrayBox& f3fab   = f3[solnLmfi];
#if (BL_SPACEDIM == 3)
        const FArrayBox& f4fab   = f4[solnLmfi];
        const FArrayBox& f5fab   = f5[solnLmfi];
#endif

#if (BL_SPACEDIM == 2)
        amrex_lp_gsrb(
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
	    tbx.loVect(), tbx.hiVect(), vbx.loVect(), vbx.hiVect(),
            &nc, h[level].data(), &redBlackFlag);
#endif

#if (BL_SPACEDIM == 3)
        amrex_lp_gsrb(
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
	    tbx.loVect(), tbx.hiVect(), vbx.loVect(), vbx.hiVect(),
	    &nc, h[level].data(), &redBlackFlag);
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
  int src_comp = 0;
  int dst_comp = 0;
  int num_comp = 1;
  Fapply(y,dst_comp,x,src_comp,num_comp,level);
}

void
Laplacian::Fapply (MultiFab&       y,
		   int             dst_comp,
                   const MultiFab& x,
		   int             src_comp,
		   int             num_comp,
                   int             level)
{
    BL_PROFILE("Laplacian::Fapply()");

    const bool tiling = true;
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter ymfi(y,tiling); ymfi.isValid(); ++ymfi)
    {
        const Box&       tbx  = ymfi.tilebox();
        FArrayBox&       yfab = y[ymfi];
        const FArrayBox& xfab = x[ymfi];

        amrex_lp_adotx(yfab.dataPtr(dst_comp), 
                   ARLIM(yfab.loVect()), ARLIM(yfab.hiVect()),
                   xfab.dataPtr(src_comp), 
                   ARLIM(xfab.loVect()), ARLIM(xfab.hiVect()),
                   tbx.loVect(), tbx.hiVect(), &num_comp,
                   h[level].data());
    }
}

}
