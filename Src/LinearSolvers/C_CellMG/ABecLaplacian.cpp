//BL_COPYRIGHT_NOTICE

//
// $Id: ABecLaplacian.cpp,v 1.9 2000-08-02 16:06:44 car Exp $
//

#include <ABecLaplacian.H>
#include <ABec_F.H>
#include <ParallelDescriptor.H>

#ifdef BL3_PROFILING
#include <BoxLib3/Profiler.H>
#endif
#ifdef BL3_PTHREADS
#include <BoxLib3/WorkQueue.H>
extern BoxLib3::WorkQueue wrkq;
#endif

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
    if (level >= a_valid.length() || a_valid[level] == false)
    {
        if (acoefs.length() < level+1)
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
    
    if (level >= b_valid.length() || b_valid[level] == false)
    {
        if (bcoefs.length() < level+1)
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
    applyBC(in,level,bc_mode);
    const BoxArray& bxa = gbox[level];
    const MultiFab& a   = aCoefficients(level);
    const MultiFab& bX  = bCoefficients(0,level);
    const MultiFab& bY  = bCoefficients(1,level);
#if (BL_SPACEDIM == 3)
    const MultiFab& bZ  = bCoefficients(2,level);
#endif
    int nc = in.nComp();

    for (MultiFabIterator inmfi(in); inmfi.isValid(); ++inmfi)
    {
        DependentMultiFabIterator amfi(inmfi,  a);
        DependentMultiFabIterator bXmfi(inmfi, bX);
        DependentMultiFabIterator bYmfi(inmfi, bY);
        DependentMultiFabIterator xflmfi(inmfi, xflux);
        DependentMultiFabIterator yflmfi(inmfi, yflux);
#if (BL_SPACEDIM == 3)	
        DependentMultiFabIterator bZmfi(inmfi, bZ);
        DependentMultiFabIterator zflmfi(inmfi, zflux);
#endif
        BL_ASSERT(bxa[inmfi.index()] == inmfi.validbox());

        FORT_FLUX(inmfi().dataPtr(),
		  ARLIM(inmfi().loVect()), ARLIM(inmfi().hiVect()),
		  &alpha, &beta, amfi().dataPtr(), 
		  ARLIM(amfi().loVect()), ARLIM(amfi().hiVect()),
		  bXmfi().dataPtr(), 
		  ARLIM(bXmfi().loVect()), ARLIM(bXmfi().hiVect()),
		  bYmfi().dataPtr(), 
		  ARLIM(bYmfi().loVect()), ARLIM(bYmfi().hiVect()),
#if (BL_SPACEDIM == 3)
		  bZmfi().dataPtr(), 
		  ARLIM(bZmfi().loVect()), ARLIM(bZmfi().hiVect()),
#endif
		  inmfi.validbox().loVect(), inmfi.validbox().hiVect(), &nc,
		  h[level],
		  xflmfi().dataPtr(),
		  ARLIM(xflmfi().loVect()), ARLIM(xflmfi().hiVect()),
		  yflmfi().dataPtr(),
		  ARLIM(yflmfi().loVect()), ARLIM(yflmfi().hiVect())
#if (BL_SPACEDIM == 3)
		  ,zflmfi().dataPtr(),
		  ARLIM(zflmfi().loVect()), ARLIM(zflmfi().hiVect())
#endif
		  );
    }
}


#ifdef BL3_PTHREADS
class task_gsrb
  : public BoxLib3::WorkQueue::task
{
public:
  task_gsrb(FArrayBox& solnL_,
	    const FArrayBox& rhsL_,
	    Real alpha_, Real beta_,
	    const FArrayBox& a_,
	    const FArrayBox& bX_,
	    const FArrayBox& bY_,
#if BL_SPACEDIM == 3
	    const FArrayBox& bZ_,
#endif
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
    : solnL(solnL_),
      rhsL(rhsL_),
      alpha(alpha_), beta(beta_),
      a(a_),
      bX(bX_),
      bY(bY_),
#if BL_SPACEDIM == 3
      bZ(bZ_),
#endif
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
  const FArrayBox& bX;
  const FArrayBox& bY;
#if BL_SPACEDIM == 3
  const FArrayBox& bZ;
#endif
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
task_gsrb::run()
{
  BL3_PROFILE(BL3_PROFILE_THIS_NAME() + "::run()");
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
#ifdef BL3_PROFILING
    BL3_PROFILE(BL3_PROFILE_THIS_NAME() + "::Fsmooth()");
#endif
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
    const MultiFab  &bX = bCoefficients(0,level);
    const MultiFab  &bY = bCoefficients(1,level);
#if (BL_SPACEDIM > 2)    
    const MultiFab  &bZ = bCoefficients(2,level);
#endif    

    int nc = solnL.nComp();
    for (MultiFabIterator solnLmfi(solnL); solnLmfi.isValid(); ++solnLmfi)
    {
        DependentMultiFabIterator rhsLmfi(solnLmfi, rhsL);
        DependentMultiFabIterator amfi(solnLmfi,  a);
        DependentMultiFabIterator bXmfi(solnLmfi, bX);
        DependentMultiFabIterator bYmfi(solnLmfi, bY);
#if (BL_SPACEDIM > 2)    
        DependentMultiFabIterator bZmfi(solnLmfi, bZ);
#endif    
        DependentFabSetIterator f0fsi(solnLmfi, f0);
        DependentFabSetIterator f1fsi(solnLmfi, f1);
        DependentFabSetIterator f2fsi(solnLmfi, f2);
        DependentFabSetIterator f3fsi(solnLmfi, f3);
#if (BL_SPACEDIM > 2)
        DependentFabSetIterator f4fsi(solnLmfi, f4);
        DependentFabSetIterator f5fsi(solnLmfi, f5);
#endif    
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

#ifdef BL3_PTHREADS
	wrkq.add(new task_gsrb(solnL[gn],
			       rhsL[gn],
			       alpha, beta,
			       a[gn],
			       bX[gn],
			       bY[gn],
#if BL_SPACEDIM == 3
			       bZ[gn],
#endif
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
        FORT_GSRB(solnLmfi().dataPtr(), ARLIM(solnLmfi().loVect()),ARLIM(solnLmfi().hiVect()),
                  rhsLmfi().dataPtr(), ARLIM(rhsLmfi().loVect()), ARLIM(rhsLmfi().hiVect()),
                  &alpha, &beta,
                  amfi().dataPtr(), ARLIM(amfi().loVect()),    ARLIM(amfi().hiVect()),
                  bXmfi().dataPtr(), ARLIM(bXmfi().loVect()),   ARLIM(bXmfi().hiVect()),
                  bYmfi().dataPtr(), ARLIM(bYmfi().loVect()),   ARLIM(bYmfi().hiVect()),
                  f0fsi().dataPtr(), ARLIM(f0fsi().loVect()),   ARLIM(f0fsi().hiVect()),
                  m0.dataPtr(), ARLIM(m0.loVect()),   ARLIM(m0.hiVect()),
                  f1fsi().dataPtr(), ARLIM(f1fsi().loVect()),   ARLIM(f1fsi().hiVect()),
                  m1.dataPtr(), ARLIM(m1.loVect()),   ARLIM(m1.hiVect()),
                  f2fsi().dataPtr(), ARLIM(f2fsi().loVect()),   ARLIM(f2fsi().hiVect()),
                  m2.dataPtr(), ARLIM(m2.loVect()),   ARLIM(m2.hiVect()),
                  f3fsi().dataPtr(), ARLIM(f3fsi().loVect()),   ARLIM(f3fsi().hiVect()),
                  m3.dataPtr(), ARLIM(m3.loVect()),   ARLIM(m3.hiVect()),
                  solnLmfi.validbox().loVect(), solnLmfi.validbox().hiVect(),
                  &nc, h[level], &redBlackFlag);
#endif

#if (BL_SPACEDIM == 3)
        FORT_GSRB(solnLmfi().dataPtr(), ARLIM(solnLmfi().loVect()),ARLIM(solnLmfi().hiVect()),
                  rhsLmfi().dataPtr(), ARLIM(rhsLmfi().loVect()), ARLIM(rhsLmfi().hiVect()),
                  &alpha, &beta,
                  amfi().dataPtr(), ARLIM(amfi().loVect()), ARLIM(amfi().hiVect()),
                  bXmfi().dataPtr(), ARLIM(bXmfi().loVect()), ARLIM(bXmfi().hiVect()),
                  bYmfi().dataPtr(), ARLIM(bYmfi().loVect()), ARLIM(bYmfi().hiVect()),
                  bZmfi().dataPtr(), ARLIM(bZmfi().loVect()), ARLIM(bZmfi().hiVect()),
                  f0fsi().dataPtr(), ARLIM(f0fsi().loVect()), ARLIM(f0fsi().hiVect()),
                  m0.dataPtr(), ARLIM(m0.loVect()), ARLIM(m0.hiVect()),
                  f1fsi().dataPtr(), ARLIM(f1fsi().loVect()), ARLIM(f1fsi().hiVect()),
                  m1.dataPtr(), ARLIM(m1.loVect()), ARLIM(m1.hiVect()),
                  f2fsi().dataPtr(), ARLIM(f2fsi().loVect()), ARLIM(f2fsi().hiVect()),
                  m2.dataPtr(), ARLIM(m2.loVect()), ARLIM(m2.hiVect()),
                  f3fsi().dataPtr(), ARLIM(f3fsi().loVect()), ARLIM(f3fsi().hiVect()),
                  m3.dataPtr(), ARLIM(m3.loVect()), ARLIM(m3.hiVect()),
                  f4fsi().dataPtr(), ARLIM(f4fsi().loVect()), ARLIM(f4fsi().hiVect()),
                  m4.dataPtr(), ARLIM(m4.loVect()), ARLIM(m4.hiVect()),
                  f5fsi().dataPtr(), ARLIM(f5fsi().loVect()), ARLIM(f5fsi().hiVect()),
                  m5.dataPtr(), ARLIM(m5.loVect()), ARLIM(m5.hiVect()),
                  solnLmfi.validbox().loVect(), solnLmfi.validbox().hiVect(),
                  &nc, h[level], &redBlackFlag);
#endif
#endif
    }
#ifdef BL3_PTHREADS
    wrkq.wait();
#endif
}

void
ABecLaplacian::Fapply (MultiFab&       y,
                       const MultiFab& x,
                       int             level)
{
#ifdef BL3_PROFILING
    BL3_PROFILE(BL3_PROFILE_THIS_NAME() + "::Fapply()");
#endif
    const BoxArray& bxa = gbox[level];
    const MultiFab& a   = aCoefficients(level);
    const MultiFab& bX  = bCoefficients(0,level);
    const MultiFab& bY  = bCoefficients(1,level);
#if (BL_SPACEDIM > 2)
    const MultiFab& bZ  = bCoefficients(2,level);
#endif
    int nc = y.nComp();

    for (MultiFabIterator ymfi(y); ymfi.isValid(); ++ymfi)
    {
        DependentMultiFabIterator xmfi(ymfi,  x);
        DependentMultiFabIterator amfi(ymfi,  a);
        DependentMultiFabIterator bXmfi(ymfi, bX);
        DependentMultiFabIterator bYmfi(ymfi, bY);
#if (BL_SPACEDIM > 2)    
        DependentMultiFabIterator bZmfi(ymfi, bZ);
#endif    
        BL_ASSERT(bxa[ymfi.index()] == ymfi.validbox());

#if (BL_SPACEDIM == 2)
        FORT_ADOTX(ymfi().dataPtr(),
                   ARLIM(ymfi().loVect()),ARLIM(ymfi().hiVect()),
                   xmfi().dataPtr(),
                   ARLIM(xmfi().loVect()), ARLIM(xmfi().hiVect()),
                   &alpha, &beta, amfi().dataPtr(), 
                   ARLIM(amfi().loVect()), ARLIM(amfi().hiVect()),
                   bXmfi().dataPtr(), 
                   ARLIM(bXmfi().loVect()), ARLIM(bXmfi().hiVect()),
                   bYmfi().dataPtr(), 
                   ARLIM(bYmfi().loVect()), ARLIM(bYmfi().hiVect()),
                   ymfi.validbox().loVect(), ymfi.validbox().hiVect(), &nc,
                   h[level]);
#endif
#if (BL_SPACEDIM ==3)
        FORT_ADOTX(ymfi().dataPtr(),
                   ARLIM(ymfi().loVect()), ARLIM(ymfi().hiVect()),
                   xmfi().dataPtr(),
                   ARLIM(xmfi().loVect()), ARLIM(xmfi().hiVect()),
                   &alpha, &beta, amfi().dataPtr(), 
                   ARLIM(amfi().loVect()), ARLIM(amfi().hiVect()),
                   bXmfi().dataPtr(), 
                   ARLIM(bXmfi().loVect()), ARLIM(bXmfi().hiVect()),
                   bYmfi().dataPtr(), 
                   ARLIM(bYmfi().loVect()), ARLIM(bYmfi().hiVect()),
                   bZmfi().dataPtr(), 
                   ARLIM(bZmfi().loVect()), ARLIM(bZmfi().hiVect()),
                   ymfi.validbox().loVect(), ymfi.validbox().hiVect(), &nc,
                   h[level]);
#endif
    }
}
