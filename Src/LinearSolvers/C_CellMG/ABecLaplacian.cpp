
//
// $Id: ABecLaplacian.cpp,v 1.17 2001-08-21 22:15:41 car Exp $
//
#include <winstd.H>

#include <algorithm>

#include <ABecLaplacian.H>
#include <ABec_F.H>
#include <ParallelDescriptor.H>
#include <Profiler.H>

#include <WorkQueue.H>

Real ABecLaplacian::a_def     = 0.0;
Real ABecLaplacian::b_def     = 1.0;
Real ABecLaplacian::alpha_def = 1.0;
Real ABecLaplacian::beta_def  = 1.0;

ABecLaplacian::ABecLaplacian (const BndryData& _bd,
                              Real             _h)
    :
    LinOp(_bd,_h),
    alpha(alpha_def),
    beta(beta_def)
{
    initCoefficients(_bd.boxes());
}

ABecLaplacian::ABecLaplacian (const BndryData& _bd,
                              const Real*      _h)
    :
    LinOp(_bd,_h),
    alpha(alpha_def),
    beta(beta_def)
{
    initCoefficients(_bd.boxes());
}

ABecLaplacian::~ABecLaplacian ()
{
    clearToLevel(-1);
}

Real
ABecLaplacian::norm(int nm, int level)
{
  BL_ASSERT(nm == 0);
  const MultiFab& a   = aCoefficients(level);

  D_TERM(const MultiFab& bX  = bCoefficients(0,level);,
         const MultiFab& bY  = bCoefficients(1,level);,
         const MultiFab& bZ  = bCoefficients(2,level););

  const int nc = a.nComp();
  Real res = 0.0;
  for (MFIter amfi(a); amfi.isValid(); ++amfi)
  {
      Real tres;
#if (BL_SPACEDIM==2)
      FORT_NORMA(&tres,
		 &alpha, &beta,
		 a[amfi].dataPtr(),  ARLIM(a[amfi].loVect()), ARLIM(a[amfi].hiVect()),
		 bX[amfi].dataPtr(), ARLIM(bX[amfi].loVect()), ARLIM(bX[amfi].hiVect()),
		 bY[amfi].dataPtr(), ARLIM(bY[amfi].loVect()), ARLIM(bY[amfi].hiVect()),
		 amfi.validbox().loVect(), amfi.validbox().hiVect(), &nc,
		 h[level]);
#elif (BL_SPACEDIM==3)

      FORT_NORMA(&tres,
		 &alpha, &beta,
		 a[amfi].dataPtr(),  ARLIM(a[amfi].loVect()), ARLIM(a[amfi].hiVect()),
		 bX[amfi].dataPtr(), ARLIM(bX[amfi].loVect()), ARLIM(bX[amfi].hiVect()),
		 bY[amfi].dataPtr(), ARLIM(bY[amfi].loVect()), ARLIM(bY[amfi].hiVect()),
		 bZ[amfi].dataPtr(), ARLIM(bZ[amfi].loVect()), ARLIM(bZ[amfi].hiVect()),
		 amfi.validbox().loVect(), amfi.validbox().hiVect(), &nc,
		 h[level]);
#endif
      res = std::max(res, tres);
  }
  ParallelDescriptor::ReduceRealMax(res);
  return res;
}

void
ABecLaplacian::clearToLevel (int level)
{
    BL_ASSERT(level >= -1);

    for (int i = level+1; i < numLevels(); ++i)
    {
        delete acoefs[i];
        a_valid[i] = false;
        for (int j = 0; j < BL_SPACEDIM; ++j)
        {
            delete bcoefs[i][j];
        }
        b_valid[i] = false;
    }
}

void
ABecLaplacian::prepareForLevel (int level)
{
    LinOp::prepareForLevel(level);

    if (level == 0 )
        return;

    prepareForLevel(level-1);
    //
    // If coefficients were marked invalid, or if not yet made, make new ones
    // (Note: makeCoefficients is a LinOp routine, and it allocates AND
    // fills coefficients.  A more efficient implementation would allocate
    // and fill in separate steps--we could then use the a_valid bool
    // along with the length of a_valid to separately determine whether to
    // fill or allocate the coefficient MultiFabs.
    //
    if (level >= a_valid.size() || a_valid[level] == false)
    {
        if (acoefs.size() < level+1)
        {
            acoefs.resize(level+1);
            acoefs[level] = new MultiFab;
        }
        else
        {
            delete acoefs[level];
            acoefs[level] = new MultiFab;
        }
        makeCoefficients(*acoefs[level], *acoefs[level-1], level);
        a_valid.resize(level+1);
        a_valid[level] = true;
    }
    
    if (level >= b_valid.size() || b_valid[level] == false)
    {
        if (bcoefs.size() < level+1)
        {
            bcoefs.resize(level+1);
            for(int i = 0; i < BL_SPACEDIM; ++i)
                bcoefs[level][i] = new MultiFab;
        }
        else
        {
            for(int i = 0; i < BL_SPACEDIM; ++i)
            {
                delete bcoefs[level][i];
                bcoefs[level][i] = new MultiFab;
            }
        }
        for (int i = 0; i < BL_SPACEDIM; ++i)
        {
            makeCoefficients(*bcoefs[level][i], *bcoefs[level-1][i], level);
        }
        b_valid.resize(level+1);
        b_valid[level] = true;
    }
}

void
ABecLaplacian::initCoefficients (const BoxArray& _ba)
{
    const int nComp=1;
    const int nGrow=0;
    acoefs.resize(1);
    bcoefs.resize(1);
    acoefs[0] = new MultiFab(_ba, nComp, nGrow, Fab_allocate);
    acoefs[0]->setVal(a_def);
    a_valid.resize(1);
    a_valid[0] = true;

    for (int i = 0; i < BL_SPACEDIM; ++i)
    {
        BoxArray edge_boxes(_ba);
        edge_boxes.surroundingNodes(i);
        bcoefs[0][i] = new MultiFab(edge_boxes, nComp, nGrow, Fab_allocate);
        bcoefs[0][i]->setVal(b_def);
    }
    b_valid.resize(1);
    b_valid[0] = true;
}

void
ABecLaplacian::invalidate_a_to_level (int lev)
{
    lev = (lev >= 0 ? lev : 0);
    for (int i = lev; i < numLevels(); i++)
        a_valid[i] = false;
}

void
ABecLaplacian::invalidate_b_to_level (int lev)
{
    lev = (lev >= 0 ? lev : 0);
    for (int i = lev; i < numLevels(); i++)
        b_valid[i] = false;
}

void
ABecLaplacian::compFlux (D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
			 MultiFab& in, const BC_Mode& bc_mode)
{
    int level = 0;
    int src_comp = 0;
    int num_comp = 1;
    applyBC(in,src_comp,num_comp,level,bc_mode);
    const BoxArray& bxa = gbox[level];
    const MultiFab& a   = aCoefficients(level);

    D_TERM(const MultiFab& bX  = bCoefficients(0,level);,
           const MultiFab& bY  = bCoefficients(1,level);,
           const MultiFab& bZ  = bCoefficients(2,level););

    int nc = in.nComp();

    for (MFIter inmfi(in); inmfi.isValid(); ++inmfi)
    {
        BL_ASSERT(bxa[inmfi.index()] == inmfi.validbox());

        FORT_FLUX(in[inmfi].dataPtr(),
		  ARLIM(in[inmfi].loVect()), ARLIM(in[inmfi].hiVect()),
		  &alpha, &beta, a[inmfi].dataPtr(), 
		  ARLIM(a[inmfi].loVect()), ARLIM(a[inmfi].hiVect()),
		  bX[inmfi].dataPtr(), 
		  ARLIM(bX[inmfi].loVect()), ARLIM(bX[inmfi].hiVect()),
		  bY[inmfi].dataPtr(), 
		  ARLIM(bY[inmfi].loVect()), ARLIM(bY[inmfi].hiVect()),
#if (BL_SPACEDIM == 3)
		  bZ[inmfi].dataPtr(), 
		  ARLIM(bZ[inmfi].loVect()), ARLIM(bZ[inmfi].hiVect()),
#endif
		  inmfi.validbox().loVect(), inmfi.validbox().hiVect(), &nc,
		  h[level],
		  xflux[inmfi].dataPtr(),
		  ARLIM(xflux[inmfi].loVect()), ARLIM(xflux[inmfi].hiVect()),
		  yflux[inmfi].dataPtr(),
		  ARLIM(yflux[inmfi].loVect()), ARLIM(yflux[inmfi].hiVect())
#if (BL_SPACEDIM == 3)
		  ,zflux[inmfi].dataPtr(),
		  ARLIM(zflux[inmfi].loVect()), ARLIM(zflux[inmfi].hiVect())
#endif
		  );
    }
}


#ifdef BL_THREADS
class task_gsrb
  :
    public WorkQueue::task
{
public:
  task_gsrb(FArrayBox& solnL_,
	    const FArrayBox& rhsL_,
	    Real alpha_, Real beta_,
	    const FArrayBox& a_,
	    D_DECL(const FArrayBox& bX_,
                   const FArrayBox& bY_,
                   const FArrayBox& bZ_),
	    const FArrayBox& f0_, const Mask& m0_,
	    const FArrayBox& f1_, const Mask& m1_,
	    const FArrayBox& f2_, const Mask& m2_,
	    const FArrayBox& f3_, const Mask& m3_,
#if BL_SPACEDIM == 3
	    const FArrayBox& f4_, const Mask& m4_,
	    const FArrayBox& f5_, const Mask& m5_,
#endif
	    const Box& vbox_,
	    int nc_,
	    const Real* h_,
	    int redBlackFlag_)
    :
      solnL(solnL_),
      rhsL(rhsL_),
      alpha(alpha_), beta(beta_),
      a(a_),
      D_DECL(bX(bX_),
             bY(bY_),
             bZ(bZ_)),
      f0(f0_), m0(m0_),
      f1(f1_), m1(m1_),
      f2(f2_), m2(m2_),
      f3(f3_), m3(m3_),
#if BL_SPACEDIM == 3
      f4(f4_), m4(m4_),
      f5(f5_), m5(m5_),
#endif
      vbox(vbox_),
      nc(nc_),
      h(h_),
      redBlackFlag(redBlackFlag_)
  {}
  virtual void run();
private:
  FArrayBox& solnL;
  const FArrayBox& rhsL;
  const Real alpha, beta;
  const FArrayBox& a;
  D_TERM(const FArrayBox& bX,
         ;const FArrayBox& bY,
         ;const FArrayBox& bZ);
  const FArrayBox& f0;
  const Mask& m0;
  const FArrayBox& f1;
  const Mask& m1;
  const FArrayBox& f2;
  const Mask& m2;
  const FArrayBox& f3;
  const Mask& m3;
#if BL_SPACEDIM == 3
  const FArrayBox& f4;
  const Mask& m4;
  const FArrayBox& f5;
  const Mask& m5;
#endif
  const Box vbox;
  const int nc;
  const Real* h;
  const int redBlackFlag;
};

void
task_gsrb::run ()
{
  BL_PROFILE(BL_PROFILE_THIS_NAME() + "::run()");

  FORT_GSRB(solnL.dataPtr(), ARLIM(solnL.loVect()),ARLIM(solnL.hiVect()),
	    rhsL.dataPtr(), ARLIM(rhsL.loVect()), ARLIM(rhsL.hiVect()),
	    &alpha, &beta,
	    a.dataPtr(), ARLIM(a.loVect()), ARLIM(a.hiVect()),
	    bX.dataPtr(), ARLIM(bX.loVect()), ARLIM(bX.hiVect()),
	    bY.dataPtr(), ARLIM(bY.loVect()), ARLIM(bY.hiVect()),
#if BL_SPACEDIM==3
	    bZ.dataPtr(), ARLIM(bZ.loVect()), ARLIM(bZ.hiVect()),
#endif
	    f0.dataPtr(), ARLIM(f0.loVect()), ARLIM(f0.hiVect()),
	    m0.dataPtr(), ARLIM(m0.loVect()), ARLIM(m0.hiVect()),
	    f1.dataPtr(), ARLIM(f1.loVect()), ARLIM(f1.hiVect()),
	    m1.dataPtr(), ARLIM(m1.loVect()), ARLIM(m1.hiVect()),
	    f2.dataPtr(), ARLIM(f2.loVect()), ARLIM(f2.hiVect()),
	    m2.dataPtr(), ARLIM(m2.loVect()), ARLIM(m2.hiVect()),
	    f3.dataPtr(), ARLIM(f3.loVect()), ARLIM(f3.hiVect()),
	    m3.dataPtr(), ARLIM(m3.loVect()), ARLIM(m3.hiVect()),
#if BL_SPACEDIM==3  
	    f4.dataPtr(), ARLIM(f4.loVect()), ARLIM(f4.hiVect()),
	    m4.dataPtr(), ARLIM(m4.loVect()), ARLIM(m4.hiVect()),
	    f5.dataPtr(), ARLIM(f5.loVect()), ARLIM(f5.hiVect()),
	    m5.dataPtr(), ARLIM(m5.loVect()), ARLIM(m5.hiVect()),
#endif
	    vbox.loVect(), vbox.hiVect(),
	    &nc, h, &redBlackFlag);
}
#endif	
        
//
// Must be defined for MultiGrid/CGSolver to work.
//

void
ABecLaplacian::Fsmooth (MultiFab&       solnL,
                        const MultiFab& rhsL,
                        int             level,
                        int             redBlackFlag)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::Fsmooth()");

    const BoxArray& bxa = gbox[level];

    OrientationIter oitr;

    const FabSet &f0 = (*undrrelxr[level])[oitr()]; oitr++;
    const FabSet &f1 = (*undrrelxr[level])[oitr()]; oitr++;
    const FabSet &f2 = (*undrrelxr[level])[oitr()]; oitr++;
    const FabSet &f3 = (*undrrelxr[level])[oitr()]; oitr++;
#if (BL_SPACEDIM > 2)
    const FabSet &f4 = (*undrrelxr[level])[oitr()]; oitr++;
    const FabSet &f5 = (*undrrelxr[level])[oitr()]; oitr++;
#endif    
    const MultiFab  &a = aCoefficients(level);
    D_TERM(const MultiFab  &bX = bCoefficients(0,level);,
           const MultiFab  &bY = bCoefficients(1,level);,
           const MultiFab  &bZ = bCoefficients(2,level););

    int nc = solnL.nComp();
    for (MFIter solnLmfi(solnL); solnLmfi.isValid(); ++solnLmfi)
    {
        oitr.rewind();
        int gn = solnLmfi.index();
        const Mask& m0 = *maskvals[level][gn][oitr()]; oitr++;
        const Mask& m1 = *maskvals[level][gn][oitr()]; oitr++;
        const Mask& m2 = *maskvals[level][gn][oitr()]; oitr++;
        const Mask& m3 = *maskvals[level][gn][oitr()]; oitr++;
#if (BL_SPACEDIM > 2)
        const Mask& m4 = *maskvals[level][gn][oitr()]; oitr++;
        const Mask& m5 = *maskvals[level][gn][oitr()]; oitr++;
#endif

        BL_ASSERT(bxa[solnLmfi.index()] == solnLmfi.validbox());

#ifdef BL_THREADS
	BoxLib::theWorkQueue().add(new task_gsrb(solnL[gn],
			       rhsL[gn],
			       alpha, beta,
			       a[gn],
			       D_DECL(bX[gn], bY[gn], bZ[gn]),
			       f0[gn], m0,
			       f1[gn], m1,
			       f2[gn], m2,
			       f3[gn], m3,
#if BL_SPACEDIM == 3
			       f4[gn], m4,
			       f5[gn], m5,
#endif
			       solnLmfi.validbox(),
			       nc, h[level], redBlackFlag));
#else
#if (BL_SPACEDIM == 2)
        FORT_GSRB(solnL[solnLmfi].dataPtr(), ARLIM(solnL[solnLmfi].loVect()),ARLIM(solnL[solnLmfi].hiVect()),
                  rhsL[solnLmfi].dataPtr(), ARLIM(rhsL[solnLmfi].loVect()), ARLIM(rhsL[solnLmfi].hiVect()),
                  &alpha, &beta,
                  a[solnLmfi].dataPtr(), ARLIM(a[solnLmfi].loVect()),    ARLIM(a[solnLmfi].hiVect()),
                  bX[solnLmfi].dataPtr(), ARLIM(bX[solnLmfi].loVect()),   ARLIM(bX[solnLmfi].hiVect()),
                  bY[solnLmfi].dataPtr(), ARLIM(bY[solnLmfi].loVect()),   ARLIM(bY[solnLmfi].hiVect()),
                  f0[solnLmfi.index()].dataPtr(), ARLIM(f0[solnLmfi.index()].loVect()),   ARLIM(f0[solnLmfi.index()].hiVect()),
                  m0.dataPtr(), ARLIM(m0.loVect()),   ARLIM(m0.hiVect()),
                  f1[solnLmfi.index()].dataPtr(), ARLIM(f1[solnLmfi.index()].loVect()),   ARLIM(f1[solnLmfi.index()].hiVect()),
                  m1.dataPtr(), ARLIM(m1.loVect()),   ARLIM(m1.hiVect()),
                  f2[solnLmfi.index()].dataPtr(), ARLIM(f2[solnLmfi.index()].loVect()),   ARLIM(f2[solnLmfi.index()].hiVect()),
                  m2.dataPtr(), ARLIM(m2.loVect()),   ARLIM(m2.hiVect()),
                  f3[solnLmfi.index()].dataPtr(), ARLIM(f3[solnLmfi.index()].loVect()),   ARLIM(f3[solnLmfi.index()].hiVect()),
                  m3.dataPtr(), ARLIM(m3.loVect()),   ARLIM(m3.hiVect()),
                  solnLmfi.validbox().loVect(), solnLmfi.validbox().hiVect(),
                  &nc, h[level], &redBlackFlag);
#endif

#if (BL_SPACEDIM == 3)
        FORT_GSRB(solnL[solnLmfi].dataPtr(), ARLIM(solnL[solnLmfi].loVect()),ARLIM(solnL[solnLmfi].hiVect()),
                  rhsL[solnLmfi].dataPtr(), ARLIM(rhsL[solnLmfi].loVect()), ARLIM(rhsL[solnLmfi].hiVect()),
                  &alpha, &beta,
                  a[solnLmfi].dataPtr(), ARLIM(a[solnLmfi].loVect()), ARLIM(a[solnLmfi].hiVect()),
                  bX[solnLmfi].dataPtr(), ARLIM(bX[solnLmfi].loVect()), ARLIM(bX[solnLmfi].hiVect()),
                  bY[solnLmfi].dataPtr(), ARLIM(bY[solnLmfi].loVect()), ARLIM(bY[solnLmfi].hiVect()),
                  bZ[solnLmfi].dataPtr(), ARLIM(bZ[solnLmfi].loVect()), ARLIM(bZ[solnLmfi].hiVect()),
                  f0[solnLmfi.index()].dataPtr(), ARLIM(f0[solnLmfi.index()].loVect()), ARLIM(f0[solnLmfi.index()].hiVect()),
                  m0.dataPtr(), ARLIM(m0.loVect()), ARLIM(m0.hiVect()),
                  f1[solnLmfi.index()].dataPtr(), ARLIM(f1[solnLmfi.index()].loVect()), ARLIM(f1[solnLmfi.index()].hiVect()),
                  m1.dataPtr(), ARLIM(m1.loVect()), ARLIM(m1.hiVect()),
                  f2[solnLmfi.index()].dataPtr(), ARLIM(f2[solnLmfi.index()].loVect()), ARLIM(f2[solnLmfi.index()].hiVect()),
                  m2.dataPtr(), ARLIM(m2.loVect()), ARLIM(m2.hiVect()),
                  f3[solnLmfi.index()].dataPtr(), ARLIM(f3[solnLmfi.index()].loVect()), ARLIM(f3[solnLmfi.index()].hiVect()),
                  m3.dataPtr(), ARLIM(m3.loVect()), ARLIM(m3.hiVect()),
                  f4[solnLmfi.index()].dataPtr(), ARLIM(f4[solnLmfi.index()].loVect()), ARLIM(f4[solnLmfi.index()].hiVect()),
                  m4.dataPtr(), ARLIM(m4.loVect()), ARLIM(m4.hiVect()),
                  f5[solnLmfi.index()].dataPtr(), ARLIM(f5[solnLmfi.index()].loVect()), ARLIM(f5[solnLmfi.index()].hiVect()),
                  m5.dataPtr(), ARLIM(m5.loVect()), ARLIM(m5.hiVect()),
                  solnLmfi.validbox().loVect(), solnLmfi.validbox().hiVect(),
                  &nc, h[level], &redBlackFlag);
#endif
#endif
    }
    BoxLib::theWorkQueue().wait();
}

void
ABecLaplacian::Fapply (MultiFab&       y,
                       const MultiFab& x,
                       int             level)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::Fapply()");

    const BoxArray& bxa = gbox[level];
    const MultiFab& a   = aCoefficients(level);
    D_TERM(const MultiFab& bX  = bCoefficients(0,level);,
           const MultiFab& bY  = bCoefficients(1,level);,
           const MultiFab& bZ  = bCoefficients(2,level););
    int nc = y.nComp();

    for (MFIter ymfi(y); ymfi.isValid(); ++ymfi)
    {
        BL_ASSERT(bxa[ymfi.index()] == ymfi.validbox());

#if (BL_SPACEDIM == 2)
        FORT_ADOTX(y[ymfi].dataPtr(),
                   ARLIM(y[ymfi].loVect()),ARLIM(y[ymfi].hiVect()),
                   x[ymfi].dataPtr(),
                   ARLIM(x[ymfi].loVect()), ARLIM(x[ymfi].hiVect()),
                   &alpha, &beta, a[ymfi].dataPtr(), 
                   ARLIM(a[ymfi].loVect()), ARLIM(a[ymfi].hiVect()),
                   bX[ymfi].dataPtr(), 
                   ARLIM(bX[ymfi].loVect()), ARLIM(bX[ymfi].hiVect()),
                   bY[ymfi].dataPtr(), 
                   ARLIM(bY[ymfi].loVect()), ARLIM(bY[ymfi].hiVect()),
                   ymfi.validbox().loVect(), ymfi.validbox().hiVect(), &nc,
                   h[level]);
#endif
#if (BL_SPACEDIM ==3)
        FORT_ADOTX(y[ymfi].dataPtr(),
                   ARLIM(y[ymfi].loVect()), ARLIM(y[ymfi].hiVect()),
                   x[ymfi].dataPtr(),
                   ARLIM(x[ymfi].loVect()), ARLIM(x[ymfi].hiVect()),
                   &alpha, &beta, a[ymfi].dataPtr(), 
                   ARLIM(a[ymfi].loVect()), ARLIM(a[ymfi].hiVect()),
                   bX[ymfi].dataPtr(), 
                   ARLIM(bX[ymfi].loVect()), ARLIM(bX[ymfi].hiVect()),
                   bY[ymfi].dataPtr(), 
                   ARLIM(bY[ymfi].loVect()), ARLIM(bY[ymfi].hiVect()),
                   bZ[ymfi].dataPtr(), 
                   ARLIM(bZ[ymfi].loVect()), ARLIM(bZ[ymfi].hiVect()),
                   ymfi.validbox().loVect(), ymfi.validbox().hiVect(), &nc,
                   h[level]);
#endif
    }
}
