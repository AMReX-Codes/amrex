
//
// $Id: DivVis.cpp,v 1.10 2001-08-01 21:51:06 lijewski Exp $
//

#include <DivVis.H>
#include <DivVis_F.H>

Real DivVis::a_def     = 0.0;
Real DivVis::b_def     = 1.0;
Real DivVis::alpha_def = 1.0;
Real DivVis::beta_def  = 1.0;

int
DivVis::numberComponents ()
{
    return BL_SPACEDIM;
}

int
DivVis::numberPhases ()
{
    return BL_SPACEDIM==2 ? 4 : 8;
}

DivVis::DivVis (const BndryData& _bd,
		Real             _h)
    :
    MCLinOp(_bd, _h),
    alpha(alpha_def),
    beta(beta_def)
{
    Real __h[BL_SPACEDIM];

    D_TERM(__h[0]=_h;, __h[1]=_h;, __h[2]=_h;);

    initConstruct(__h);
}

DivVis::DivVis (const BndryData& _bd,
	        const Real*      _h)
    :
    MCLinOp(_bd, _h),
    alpha(alpha_def),
    beta(beta_def)
{
    initConstruct(_h);
}

void
DivVis::initConstruct (const Real* _h)
{
    const int level       = 0;
    const BoxArray& grids = gbox[level];

    initCoefficients(gbox[level]);

    numcomp  = numberComponents(); // wyc
    numphase = numberPhases();     // wyc

    undrrelxr.resize(1);
    undrrelxr[level] = new BndryRegister(gbox[level], 1, 0, 0, numcomp);
    tangderiv.resize(1);
#if BL_SPACEDIM==2
    tangderiv[level] = new BndryRegister(gbox[level], 0, 1, 0, numcomp);
#elif BL_SPACEDIM==3
    tangderiv[level] = new BndryRegister(gbox[level], 0, 1, 0, numcomp*(1+3));
#else
# error "BL_SPACEDIME must be 2 or 3"
#endif
}

DivVis::~DivVis ()
{
    clearToLevel(-1);
}

void
DivVis::setScalars (Real _alpha,
		    Real _beta)
{
    alpha = _alpha;
    beta  = _beta;
}

void
DivVis::clearToLevel (int level)
{
    BL_ASSERT(level >= -1);

    for (int i = level+1; i < numLevels(); ++i)
    {
	delete acoefs[i];
	for (int j = 0; j < BL_SPACEDIM; ++j)
        {
	    delete bcoefs[i][j];
	}
    }
}

void
DivVis::prepareForLevel (int level)
{
    MCLinOp::prepareForLevel(level);
    if (level == 0)
        return;
    prepareForLevel(level-1);
    //
    // If coefficients were marked invalid, or if not yet made, make new ones
    // (Note: makeCoefficients is a MCLinOp routine, and it allocates AND
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
	    for (int i = 0; i < BL_SPACEDIM; ++i)
		bcoefs[level][i] = new MultiFab;
	}
        else
        {
	    for (int i = 0; i < BL_SPACEDIM; ++i)
            {
		delete bcoefs[level][i];
		bcoefs[level][i] = new MultiFab;
	    }
	}
	for (int i = 0; i < BL_SPACEDIM; ++i)
	    makeCoefficients(*bcoefs[level][i], *bcoefs[level-1][i], level);

	b_valid.resize(level+1);
	b_valid[level] = true;
    }
}

void
DivVis::initCoefficients (const BoxArray &_ba)
{
    const int nGrow = 0;
    const int level = 0;

    acoefs.resize(1);
    bcoefs.resize(1);
    //
    // In 2D, need 2 components for "a" to handle r-z properly (will need three
    // for r-theta-phi, but allowing only 3D cartesian for now).
    //
    const int nComp = (BL_SPACEDIM == 2  ?  2  :  1);

#ifndef NDEBUG
    if (BL_SPACEDIM == 3)
	BL_ASSERT(geomarray[level].IsCartesian());
#endif

    acoefs[level] = new MultiFab(_ba, nComp, nGrow);
    acoefs[level]->setVal(a_def);
    a_valid.resize(1);
    a_valid[level] = true;

    for (int i = 0; i < BL_SPACEDIM; ++i)
    {
	BoxArray edge_boxes(_ba);
	edge_boxes.surroundingNodes(i);
	bcoefs[level][i] = new MultiFab(edge_boxes, 1, nGrow);
	bcoefs[level][i]->setVal(b_def);
    }
    b_valid.resize(1);
    b_valid[level] = true;
}

void
DivVis::invalidate_a_to_level (int lev)
{
    lev = (lev >= 0 ? lev : 0);
    for (int i = lev; i < numLevels(); i++)
        a_valid[i]=false;
}

void
DivVis::invalidate_b_to_level (int lev)
{
    lev = (lev >= 0 ? lev : 0);
    for (int i = lev; i < numLevels(); i++)
        b_valid[i]=false;
}


void
DivVis::aCoefficients (const MultiFab& _a)
{
    BL_ASSERT(_a.ok());
    BL_ASSERT(_a.boxArray() == (acoefs[0])->boxArray());
    const int nComp = (BL_SPACEDIM == 2  ?  2  :  1);
    BL_ASSERT(_a.nComp() == nComp);
    invalidate_a_to_level(0);
    (*acoefs[0]).copy(_a,0,0,nComp);
}

void
DivVis::bCoefficients (const MultiFab& _b,
                       int             dir)
{
    BL_ASSERT(_b.ok());
    BL_ASSERT(_b.boxArray() == (bcoefs[0][dir])->boxArray());
    BL_ASSERT(_b.nComp() == 1);
    invalidate_b_to_level(0);
    (*bcoefs[0][dir]).copy(_b,0,0,1);
}

const MultiFab&
DivVis::aCoefficients (int level)
{
    prepareForLevel(level);
    return *acoefs[level];
}

const MultiFab&
DivVis::bCoefficients (int dir,
                       int level)
{
    prepareForLevel(level);
    return *bcoefs[level][dir];
}

//
// Must be defined for MultiGrid/CGSolver to work.
//
void
DivVis::Fsmooth (MultiFab&       solnL,
                 const MultiFab& rhsL,
                 int             level,
                 int             phaseflag)
{
    OrientationIter oitr;

    const FabSet& fw  = (*undrrelxr[level])[oitr()]; 
    const FabSet& tdw = (*tangderiv[level])[oitr()];
    oitr++;
    const FabSet& fs  = (*undrrelxr[level])[oitr()]; 
    const FabSet& tds = (*tangderiv[level])[oitr()];
    oitr++;
#if BL_SPACEDIM>2
    const FabSet& fb  = (*undrrelxr[level])[oitr()]; 
    const FabSet& tdb = (*tangderiv[level])[oitr()];
    oitr++;
#endif
    const FabSet& fe  = (*undrrelxr[level])[oitr()]; 
    const FabSet& tde = (*tangderiv[level])[oitr()];
    oitr++;
    const FabSet& fn  = (*undrrelxr[level])[oitr()]; 
    const FabSet& tdn = (*tangderiv[level])[oitr()];
    oitr++;
#if BL_SPACEDIM>2
    const FabSet& ft  = (*undrrelxr[level])[oitr()]; 
    const FabSet& tdt = (*tangderiv[level])[oitr()];
    oitr++;
#endif
    const MultiFab& a  = aCoefficients(level);

    D_TERM(const MultiFab& bX = bCoefficients(0,level);,
           const MultiFab& bY = bCoefficients(1,level);,
           const MultiFab& bZ = bCoefficients(2,level););

    int nc = solnL.nComp();

    for (MFIter solnLmfi(solnL); solnLmfi.isValid(); ++solnLmfi)
    {
	oitr.rewind();
        const int gn = solnLmfi.index();

	D_TERM(const Mask& mw = *maskvals[level][gn][oitr()]; oitr++;,
               const Mask& ms = *maskvals[level][gn][oitr()]; oitr++;,
               const Mask& mb = *maskvals[level][gn][oitr()]; oitr++;);

	D_TERM(const Mask& me = *maskvals[level][gn][oitr()]; oitr++;,
               const Mask& mn = *maskvals[level][gn][oitr()]; oitr++;,
               const Mask& mt = *maskvals[level][gn][oitr()]; oitr++;);

	FORT_GSRB(
	    solnL[solnLmfi].dataPtr(), 
            ARLIM(solnL[solnLmfi].loVect()),ARLIM(solnL[solnLmfi].hiVect()),
	    rhsL[solnLmfi].dataPtr(),
            ARLIM(rhsL[solnLmfi].loVect()), ARLIM(rhsL[solnLmfi].hiVect()),
	    &alpha, &beta,
	    a[solnLmfi].dataPtr(),
            ARLIM(a[solnLmfi].loVect()),    ARLIM(a[solnLmfi].hiVect()),
	    bX[solnLmfi].dataPtr(),
            ARLIM(bX[solnLmfi].loVect()),   ARLIM(bX[solnLmfi].hiVect()),
	    bY[solnLmfi].dataPtr(),
            ARLIM(bY[solnLmfi].loVect()),   ARLIM(bY[solnLmfi].hiVect()),
#if BL_SPACEDIM>2
	    bZ[solnLmfi].dataPtr(),
            ARLIM(bZ[solnLmfi].loVect()),   ARLIM(bZ[solnLmfi].hiVect()),
#endif
	    mn.dataPtr(),
	    ARLIM(mn.loVect()),ARLIM(mn.hiVect()),
	    fn[solnLmfi].dataPtr(),
            ARLIM(fn[solnLmfi].loVect()),   ARLIM(fn[solnLmfi].hiVect()),
	    me.dataPtr(),
	    ARLIM(me.loVect()),ARLIM(me.hiVect()),
	    fe[solnLmfi].dataPtr(),
            ARLIM(fe[solnLmfi].loVect()),   ARLIM(fe[solnLmfi].hiVect()),
	    mw.dataPtr(),
	    ARLIM(mw.loVect()),ARLIM(mw.hiVect()),
	    fw[solnLmfi].dataPtr(),
            ARLIM(fw[solnLmfi].loVect()),   ARLIM(fw[solnLmfi].hiVect()),
	    ms.dataPtr(),
	    ARLIM(ms.loVect()),ARLIM(ms.hiVect()),
	    fs[solnLmfi].dataPtr(),
            ARLIM(fs[solnLmfi].loVect()),   ARLIM(fs[solnLmfi].hiVect()),
#if BL_SPACEDIM>2
	    mt.dataPtr(),
	    ARLIM(mt.loVect()),ARLIM(mt.hiVect()),
	    ft[solnLmfi].dataPtr(),
            ARLIM(ft[solnLmfi].loVect()),   ARLIM(ft[solnLmfi].hiVect()),
	    mb.dataPtr(),
	    ARLIM(mb.loVect()),ARLIM(mb.hiVect()),
	    fb[solnLmfi].dataPtr(),
            ARLIM(fb[solnLmfi].loVect()),   ARLIM(fb[solnLmfi].hiVect()),
#endif
	    tdn[solnLmfi].dataPtr(),
	    ARLIM(tdn[solnLmfi].loVect()),ARLIM(tdn[solnLmfi].hiVect()),
	    tde[solnLmfi].dataPtr(),
	    ARLIM(tde[solnLmfi].loVect()),ARLIM(tde[solnLmfi].hiVect()),
	    tdw[solnLmfi].dataPtr(),
	    ARLIM(tdw[solnLmfi].loVect()),ARLIM(tdw[solnLmfi].hiVect()),
	    tds[solnLmfi].dataPtr(),
	    ARLIM(tds[solnLmfi].loVect()),ARLIM(tds[solnLmfi].hiVect()),
#if BL_SPACEDIM>2
	    tdt[solnLmfi].dataPtr(),
	    ARLIM(tdt[solnLmfi].loVect()),ARLIM(tdt[solnLmfi].hiVect()),
	    tdb[solnLmfi].dataPtr(),
	    ARLIM(tdb[solnLmfi].loVect()),ARLIM(tdb[solnLmfi].hiVect()),
#endif
	    solnLmfi.validbox().loVect(), solnLmfi.validbox().hiVect(),
	    h[level], nc, phaseflag);
    }
}

void 
DivVis::compFlux (D_DECL(MultiFab& xflux, 
                         MultiFab& yflux, 
                         MultiFab& zflux), 
		  MultiFab& x)
{
    const int level   = 0;
    MCBC_Mode bc_mode = MCInhomogeneous_BC;
    applyBC(x,level,bc_mode);
    
    const MultiFab& a  = aCoefficients(level);

    D_TERM(const MultiFab& bX = bCoefficients(0,level);,
           const MultiFab& bY = bCoefficients(1,level);,
           const MultiFab& bZ = bCoefficients(2,level););

    OrientationIter oitr;

    D_TERM(const FabSet& tdw = (*tangderiv[level])[oitr()]; oitr++;,
           const FabSet& tds = (*tangderiv[level])[oitr()]; oitr++;,
           const FabSet& tdb = (*tangderiv[level])[oitr()]; oitr++;);

    D_TERM(const FabSet& tde = (*tangderiv[level])[oitr()]; oitr++;,
           const FabSet& tdn = (*tangderiv[level])[oitr()]; oitr++;,
           const FabSet& tdt = (*tangderiv[level])[oitr()]; oitr++;);

    const int nc = x.nComp();

    BL_ASSERT(nc == BL_SPACEDIM);
    BL_ASSERT(nc == xflux.nComp());
    BL_ASSERT(nc == yflux.nComp());

    for (MFIter xmfi(x); xmfi.isValid(); ++xmfi)
    {
	oitr.rewind();
        const int gn = xmfi.index();

	D_TERM(const Mask& mw = *maskvals[level][gn][oitr()]; oitr++;,
               const Mask& ms = *maskvals[level][gn][oitr()]; oitr++;,
               const Mask& mb = *maskvals[level][gn][oitr()]; oitr++;);

	D_TERM(const Mask& me = *maskvals[level][gn][oitr()]; oitr++;,
	       const Mask& mn = *maskvals[level][gn][oitr()]; oitr++;,
               const Mask& mt = *maskvals[level][gn][oitr()]; oitr++;);

	FORT_DVFLUX(
	    x[xmfi].dataPtr(), 
	    ARLIM(x[xmfi].loVect()), ARLIM(x[xmfi].hiVect()),
	    &alpha, &beta,
	    a[xmfi].dataPtr(), 
	    ARLIM(a[xmfi].loVect()), ARLIM(a[xmfi].hiVect()),
	    bX[xmfi].dataPtr(), 
	    ARLIM(bX[xmfi].loVect()), ARLIM(bX[xmfi].hiVect()),
	    bY[xmfi].dataPtr(), 
	    ARLIM(bY[xmfi].loVect()), ARLIM(bY[xmfi].hiVect()),
#if BL_SPACEDIM>2
	    bZ[xmfi].dataPtr(), 
	    ARLIM(bZ[xmfi].loVect()), ARLIM(bZ[xmfi].hiVect()),
#endif
	    xflux[xmfi].dataPtr(), 
	    ARLIM(xflux[xmfi].loVect()), ARLIM(xflux[xmfi].hiVect()),
	    yflux[xmfi].dataPtr(), 
	    ARLIM(yflux[xmfi].loVect()), ARLIM(yflux[xmfi].hiVect()),
#if BL_SPACEDIM>2
	    zflux[xmfi].dataPtr(), 
	    ARLIM(zflux[xmfi].loVect()), ARLIM(zflux[xmfi].hiVect()),
#endif
	    mn.dataPtr(),
	    ARLIM(mn.loVect()),ARLIM(mn.hiVect()),
	    me.dataPtr(),
	    ARLIM(me.loVect()),ARLIM(me.hiVect()),
	    mw.dataPtr(),
	    ARLIM(mw.loVect()),ARLIM(mw.hiVect()),
	    ms.dataPtr(),
	    ARLIM(ms.loVect()),ARLIM(ms.hiVect()),
#if BL_SPACEDIM>2
	    mt.dataPtr(),
	    ARLIM(mt.loVect()),ARLIM(mt.hiVect()),
	    mb.dataPtr(),
	    ARLIM(mb.loVect()),ARLIM(mb.hiVect()),
#endif
	    tdn[xmfi].dataPtr(),
	    ARLIM(tdn[xmfi].loVect()),ARLIM(tdn[xmfi].hiVect()),
	    tde[xmfi].dataPtr(),
	    ARLIM(tde[xmfi].loVect()),ARLIM(tde[xmfi].hiVect()),
	    tdw[xmfi].dataPtr(),
	    ARLIM(tdw[xmfi].loVect()),ARLIM(tdw[xmfi].hiVect()),
	    tds[xmfi].dataPtr(),
	    ARLIM(tds[xmfi].loVect()),ARLIM(tds[xmfi].hiVect()),
#if BL_SPACEDIM>2
	    tdt[xmfi].dataPtr(),
	    ARLIM(tdt[xmfi].loVect()),ARLIM(tdt[xmfi].hiVect()),
	    tdb[xmfi].dataPtr(),
	    ARLIM(tdb[xmfi].loVect()),ARLIM(tdb[xmfi].hiVect()),
#endif
	    xmfi.validbox().loVect(), xmfi.validbox().hiVect(),
	    h[level]);
    }
}

void
DivVis::Fapply (MultiFab&       y,
                const MultiFab& x,
                int             level)
{
    const MultiFab& a = aCoefficients(level);

    D_TERM(const MultiFab& bX = bCoefficients(0,level);,
           const MultiFab& bY = bCoefficients(1,level);,
           const MultiFab& bZ = bCoefficients(2,level););

    OrientationIter oitr;

    D_TERM(const FabSet& tdw = (*tangderiv[level])[oitr()]; oitr++;,
           const FabSet& tds = (*tangderiv[level])[oitr()]; oitr++;,
           const FabSet& tdb = (*tangderiv[level])[oitr()]; oitr++;);

    D_TERM(const FabSet& tde = (*tangderiv[level])[oitr()]; oitr++;,
           const FabSet& tdn = (*tangderiv[level])[oitr()]; oitr++;,
           const FabSet& tdt = (*tangderiv[level])[oitr()]; oitr++;);

    const int nc = y.nComp();

    for (MFIter xmfi(x); xmfi.isValid(); ++xmfi)
    {
        oitr.rewind();
        const int gn = xmfi.index();

	D_TERM(const Mask& mw = *maskvals[level][gn][oitr()]; oitr++;,
               const Mask& ms = *maskvals[level][gn][oitr()]; oitr++;,
               const Mask& mb = *maskvals[level][gn][oitr()]; oitr++;);

	D_TERM(const Mask& me = *maskvals[level][gn][oitr()]; oitr++;,
               const Mask& mn = *maskvals[level][gn][oitr()]; oitr++;,
               const Mask& mt = *maskvals[level][gn][oitr()]; oitr++;);

	FORT_DVAPPLY(
	    x[xmfi].dataPtr(), 
            ARLIM(x[xmfi].loVect()), ARLIM(x[xmfi].hiVect()),
	    &alpha, &beta,
	    a[xmfi].dataPtr(), 
            ARLIM(a[xmfi].loVect()), ARLIM(a[xmfi].hiVect()),
	    bX[xmfi].dataPtr(), 
            ARLIM(bX[xmfi].loVect()), ARLIM(bX[xmfi].hiVect()),
	    bY[xmfi].dataPtr(), 
            ARLIM(bY[xmfi].loVect()), ARLIM(bY[xmfi].hiVect()),
#if BL_SPACEDIM>2
	    bZ[xmfi].dataPtr(), 
            ARLIM(bZ[xmfi].loVect()), ARLIM(bZ[xmfi].hiVect()),
#endif
	    y[xmfi].dataPtr(), 
            ARLIM(y[xmfi].loVect()), ARLIM(y[xmfi].hiVect()),
	    mn.dataPtr(),
	    ARLIM(mn.loVect()),ARLIM(mn.hiVect()),
	    me.dataPtr(),
	    ARLIM(me.loVect()),ARLIM(me.hiVect()),
	    mw.dataPtr(),
	    ARLIM(mw.loVect()),ARLIM(mw.hiVect()),
	    ms.dataPtr(),
	    ARLIM(ms.loVect()),ARLIM(ms.hiVect()),
#if BL_SPACEDIM>2
	    mt.dataPtr(),
	    ARLIM(mt.loVect()),ARLIM(mt.hiVect()),
	    mb.dataPtr(),
	    ARLIM(mb.loVect()),ARLIM(mb.hiVect()),
#endif
	    tdn[xmfi].dataPtr(),
	    ARLIM(tdn[xmfi].loVect()),ARLIM(tdn[xmfi].hiVect()),
	    tde[xmfi].dataPtr(),
	    ARLIM(tde[xmfi].loVect()),ARLIM(tde[xmfi].hiVect()),
	    tdw[xmfi].dataPtr(),
	    ARLIM(tdw[xmfi].loVect()),ARLIM(tdw[xmfi].hiVect()),
	    tds[xmfi].dataPtr(),
	    ARLIM(tds[xmfi].loVect()),ARLIM(tds[xmfi].hiVect()),
#if BL_SPACEDIM>2
	    tdt[xmfi].dataPtr(),
	    ARLIM(tdt[xmfi].loVect()),ARLIM(tdt[xmfi].hiVect()),
	    tdb[xmfi].dataPtr(),
	    ARLIM(tdb[xmfi].loVect()),ARLIM(tdb[xmfi].hiVect()),
#endif
            xmfi.validbox().loVect(), xmfi.validbox().hiVect(),
	    h[level]);
    }
}
