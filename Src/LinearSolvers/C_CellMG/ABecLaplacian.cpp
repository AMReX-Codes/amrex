//BL_COPYRIGHT_NOTICE

//
// $Id: ABecLaplacian.cpp,v 1.5 1999-01-04 18:09:00 marc Exp $
//

#include <ABecLaplacian.H>
#include <ABec_F.H>
#include <ParallelDescriptor.H>

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
    assert(level >= -1);

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
        assert(bxa[inmfi.index()] == inmfi.validbox());

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


//
// Must be defined for MultiGrid/CGSolver to work.
//

void
ABecLaplacian::Fsmooth (MultiFab&       solnL,
                        const MultiFab& rhsL,
                        int             level,
                        int             redBlackFlag)
{
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
    for (MultiFabIterator solnLmfi(solnL); solnLmfi.isValid();
         ++solnLmfi)
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

        assert(bxa[solnLmfi.index()] == solnLmfi.validbox());
        
#if (BL_SPACEDIM == 2)
        FORT_GSRB(solnLmfi().dataPtr(), 
                  ARLIM(solnLmfi().loVect()),ARLIM(solnLmfi().hiVect()),
                  rhsLmfi().dataPtr(), 
                  ARLIM(rhsLmfi().loVect()), ARLIM(rhsLmfi().hiVect()),
                  &alpha, &beta,
                  amfi().dataPtr(),  
                  ARLIM(amfi().loVect()),    ARLIM(amfi().hiVect()),
                  bXmfi().dataPtr(), 
                  ARLIM(bXmfi().loVect()),   ARLIM(bXmfi().hiVect()),
                  bYmfi().dataPtr(),
                  ARLIM(bYmfi().loVect()),   ARLIM(bYmfi().hiVect()),
                  f0fsi().dataPtr(), 
                  ARLIM(f0fsi().loVect()),   ARLIM(f0fsi().hiVect()),
                  m0.dataPtr(), 
                  ARLIM(m0.loVect()),   ARLIM(m0.hiVect()),
                  f1fsi().dataPtr(), 
                  ARLIM(f1fsi().loVect()),   ARLIM(f1fsi().hiVect()),
                  m1.dataPtr(), 
                  ARLIM(m1.loVect()),   ARLIM(m1.hiVect()),
                  f2fsi().dataPtr(), 
                  ARLIM(f2fsi().loVect()),   ARLIM(f2fsi().hiVect()),
                  m2.dataPtr(), 
                  ARLIM(m2.loVect()),   ARLIM(m2.hiVect()),
                  f3fsi().dataPtr(), 
                  ARLIM(f3fsi().loVect()),   ARLIM(f3fsi().hiVect()),
                  m3.dataPtr(), 
                  ARLIM(m3.loVect()),   ARLIM(m3.hiVect()),
                  solnLmfi.validbox().loVect(), solnLmfi.validbox().hiVect(),
                  &nc, h[level], &redBlackFlag);
#endif

#if (BL_SPACEDIM == 3)
        FORT_GSRB(solnLmfi().dataPtr(), 
                  ARLIM(solnLmfi().loVect()),ARLIM(solnLmfi().hiVect()),
                  rhsLmfi().dataPtr(), 
                  ARLIM(rhsLmfi().loVect()), ARLIM(rhsLmfi().hiVect()),
                  &alpha, &beta,
                  amfi().dataPtr(), 
                  ARLIM(amfi().loVect()), ARLIM(amfi().hiVect()),
                  bXmfi().dataPtr(), 
                  ARLIM(bXmfi().loVect()), ARLIM(bXmfi().hiVect()),
                  bYmfi().dataPtr(), 
                  ARLIM(bYmfi().loVect()), ARLIM(bYmfi().hiVect()),
                  bZmfi().dataPtr(),
                  ARLIM(bZmfi().loVect()), ARLIM(bZmfi().hiVect()),
                  f0fsi().dataPtr(), 
                  ARLIM(f0fsi().loVect()), ARLIM(f0fsi().hiVect()),
                  m0.dataPtr(), 
                  ARLIM(m0.loVect()), ARLIM(m0.hiVect()),
                  f1fsi().dataPtr(), 
                  ARLIM(f1fsi().loVect()), ARLIM(f1fsi().hiVect()),
                  m1.dataPtr(), 
                  ARLIM(m1.loVect()), ARLIM(m1.hiVect()),
                  f2fsi().dataPtr(), 
                  ARLIM(f2fsi().loVect()), ARLIM(f2fsi().hiVect()),
                  m2.dataPtr(), 
                  ARLIM(m2.loVect()), ARLIM(m2.hiVect()),
                  f3fsi().dataPtr(), 
                  ARLIM(f3fsi().loVect()), ARLIM(f3fsi().hiVect()),
                  m3.dataPtr(), 
                  ARLIM(m3.loVect()), ARLIM(m3.hiVect()),
                  f4fsi().dataPtr(), 
                  ARLIM(f4fsi().loVect()), ARLIM(f4fsi().hiVect()),
                  m4.dataPtr(), 
                  ARLIM(m4.loVect()), ARLIM(m4.hiVect()),
                  f5fsi().dataPtr(), 
                  ARLIM(f5fsi().loVect()), ARLIM(f5fsi().hiVect()),
                  m5.dataPtr(), 
                  ARLIM(m5.loVect()), ARLIM(m5.hiVect()),
                  solnLmfi.validbox().loVect(), solnLmfi.validbox().hiVect(),
                  &nc, h[level], &redBlackFlag);
#endif
    }
}

void
ABecLaplacian::Fapply (MultiFab&       y,
                       const MultiFab& x,
                       int             level)
{
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
        assert(bxa[ymfi.index()] == ymfi.validbox());

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
