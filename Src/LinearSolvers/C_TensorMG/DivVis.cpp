
//
// $Id: DivVis.cpp,v 1.8 2000-10-02 20:53:39 lijewski Exp $
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
#if BL_SPACEDIM==2
    return 4;
#else
    return 8;
#endif
}

DivVis::DivVis (const BndryData& _bd,
		Real             _h)
    :
    MCLinOp(_bd, _h),
    alpha(alpha_def),
    beta(beta_def)
{
    Real __h[BL_SPACEDIM];

    D_TERM(__h[0] = _h;,
           __h[1] = _h;,
           __h[2] = _h;);

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

    acoefs[level] = new MultiFab(_ba, nComp, nGrow, Fab_allocate);
    acoefs[level]->setVal(a_def);
    a_valid.resize(1);
    a_valid[level] = true;

    for (int i = 0; i < BL_SPACEDIM; ++i)
    {
	BoxArray edge_boxes(_ba);
	edge_boxes.surroundingNodes(i);
	bcoefs[level][i] = new MultiFab(edge_boxes, 1, nGrow, Fab_allocate);
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
    const MultiFab& bX = bCoefficients(0,level);
    const MultiFab& bY = bCoefficients(1,level);
#if BL_SPACEDIM>2
    const MultiFab  &bZ = bCoefficients(2,level);
#endif
    int nc = solnL.nComp();

    for (MultiFabIterator solnLmfi(solnL); solnLmfi.isValid(); ++solnLmfi)
    {
	DependentMultiFabIterator rhsLmfi(solnLmfi, rhsL);
	DependentMultiFabIterator amfi(solnLmfi,  a);
	DependentMultiFabIterator bXmfi(solnLmfi, bX);
	DependentMultiFabIterator bYmfi(solnLmfi, bY);
#if BL_SPACEDIM > 2
	DependentMultiFabIterator bZmfi(solnLmfi, bZ);
#endif    
	DependentFabSetIterator fwfsi(solnLmfi,  fw);
	DependentFabSetIterator tdwfsi(solnLmfi, tdw);
	DependentFabSetIterator fsfsi(solnLmfi,  fs);
	DependentFabSetIterator tdsfsi(solnLmfi, tds);
#if BL_SPACEDIM > 2
	DependentFabSetIterator fbfsi(solnLmfi,  fb);
	DependentFabSetIterator tdbfsi(solnLmfi, tdb);
#endif
	DependentFabSetIterator fefsi(solnLmfi,  fe);
	DependentFabSetIterator tdefsi(solnLmfi, tde);
	DependentFabSetIterator fnfsi(solnLmfi,  fn);
	DependentFabSetIterator tdnfsi(solnLmfi, tdn);
#if BL_SPACEDIM > 2
	DependentFabSetIterator ftfsi(solnLmfi,  ft);
	DependentFabSetIterator tdtfsi(solnLmfi, tdt);
#endif

	oitr.rewind();
        const int gn = solnLmfi.index();

	const Mask& mw = *maskvals[level][gn][oitr()]; oitr++;
	const Mask& ms = *maskvals[level][gn][oitr()]; oitr++;
#if BL_SPACEDIM>2
	const Mask& mb = *maskvals[level][gn][oitr()]; oitr++;
#endif
	const Mask& me = *maskvals[level][gn][oitr()]; oitr++;
	const Mask& mn = *maskvals[level][gn][oitr()]; oitr++;
#if BL_SPACEDIM>2
	const Mask& mt = *maskvals[level][gn][oitr()]; oitr++;
#endif

	FORT_GSRB(
	    solnLmfi().dataPtr(), 
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
#if BL_SPACEDIM>2
	    bZmfi().dataPtr(),
            ARLIM(bZmfi().loVect()),   ARLIM(bZmfi().hiVect()),
#endif
	    mn.dataPtr(),
	    ARLIM(mn.loVect()),ARLIM(mn.hiVect()),
	    fnfsi().dataPtr(),
            ARLIM(fnfsi().loVect()),   ARLIM(fnfsi().hiVect()),
	    me.dataPtr(),
	    ARLIM(me.loVect()),ARLIM(me.hiVect()),
	    fefsi().dataPtr(),
            ARLIM(fefsi().loVect()),   ARLIM(fefsi().hiVect()),
	    mw.dataPtr(),
	    ARLIM(mw.loVect()),ARLIM(mw.hiVect()),
	    fwfsi().dataPtr(),
            ARLIM(fwfsi().loVect()),   ARLIM(fwfsi().hiVect()),
	    ms.dataPtr(),
	    ARLIM(ms.loVect()),ARLIM(ms.hiVect()),
	    fsfsi().dataPtr(),
            ARLIM(fsfsi().loVect()),   ARLIM(fsfsi().hiVect()),
#if BL_SPACEDIM>2
	    mt.dataPtr(),
	    ARLIM(mt.loVect()),ARLIM(mt.hiVect()),
	    ftfsi().dataPtr(),
            ARLIM(ftfsi().loVect()),   ARLIM(ftfsi().hiVect()),
	    mb.dataPtr(),
	    ARLIM(mb.loVect()),ARLIM(mb.hiVect()),
	    fbfsi().dataPtr(),
            ARLIM(fbfsi().loVect()),   ARLIM(fbfsi().hiVect()),
#endif
	    tdnfsi().dataPtr(),
	    ARLIM(tdnfsi().loVect()),ARLIM(tdnfsi().hiVect()),
	    tdefsi().dataPtr(),
	    ARLIM(tdefsi().loVect()),ARLIM(tdefsi().hiVect()),
	    tdwfsi().dataPtr(),
	    ARLIM(tdwfsi().loVect()),ARLIM(tdwfsi().hiVect()),
	    tdsfsi().dataPtr(),
	    ARLIM(tdsfsi().loVect()),ARLIM(tdsfsi().hiVect()),
#if BL_SPACEDIM>2
	    tdtfsi().dataPtr(),
	    ARLIM(tdtfsi().loVect()),ARLIM(tdtfsi().hiVect()),
	    tdbfsi().dataPtr(),
	    ARLIM(tdbfsi().loVect()),ARLIM(tdbfsi().hiVect()),
#endif
	    solnLmfi.validbox().loVect(), solnLmfi.validbox().hiVect(),
	    h[level], nc, phaseflag);
    }
}

void 
DivVis::compFlux (MultiFab& xflux, 
		  MultiFab& yflux, 
#if BL_SPACEDIM>2
		  MultiFab& zflux, 
#endif
		  MultiFab& x)
{
    const int level   = 0;
    MCBC_Mode bc_mode = MCInhomogeneous_BC;
    applyBC(x,level,bc_mode);
    
    const MultiFab& a  = aCoefficients(level);
    const MultiFab& bX = bCoefficients(0,level);
    const MultiFab& bY = bCoefficients(1,level);
#if BL_SPACEDIM>2
    const MultiFab& bZ = bCoefficients(2,level);
#endif
    OrientationIter oitr;
    const FabSet& tdw = (*tangderiv[level])[oitr()]; oitr++;
    const FabSet& tds = (*tangderiv[level])[oitr()]; oitr++;
#if BL_SPACEDIM>2
    const FabSet& tdb = (*tangderiv[level])[oitr()]; oitr++;
#endif
    const FabSet& tde = (*tangderiv[level])[oitr()]; oitr++;
    const FabSet& tdn = (*tangderiv[level])[oitr()]; oitr++;
#if BL_SPACEDIM>2
    const FabSet& tdt = (*tangderiv[level])[oitr()]; oitr++;
#endif
    const int nc = x.nComp();

    BL_ASSERT(nc == BL_SPACEDIM);
    BL_ASSERT(nc == xflux.nComp());
    BL_ASSERT(nc == yflux.nComp());

    for (MultiFabIterator xmfi(x); xmfi.isValid(); ++xmfi)
    {
	DependentMultiFabIterator amfi(xmfi,  a);
	DependentMultiFabIterator bXmfi(xmfi, bX);
	DependentMultiFabIterator xfluxmfi(xmfi, xflux);
	DependentMultiFabIterator bYmfi(xmfi, bY);
	DependentMultiFabIterator yfluxmfi(xmfi, yflux);
#if BL_SPACEDIM > 2
	DependentMultiFabIterator bZmfi(xmfi, bZ);
	DependentMultiFabIterator zfluxmfi(xmfi, zflux);
#endif    
	DependentFabSetIterator tdwfsi(xmfi, tdw);
	DependentFabSetIterator tdsfsi(xmfi, tds);
#if BL_SPACEDIM > 2
	DependentFabSetIterator tdbfsi(xmfi, tdb);
#endif
	DependentFabSetIterator tdefsi(xmfi, tde);
	DependentFabSetIterator tdnfsi(xmfi, tdn);
#if BL_SPACEDIM > 2
	DependentFabSetIterator tdtfsi(xmfi, tdt);
#endif

	oitr.rewind();
        const int gn = xmfi.index();
	const Mask& mw = *maskvals[level][gn][oitr()]; oitr++;
	const Mask& ms = *maskvals[level][gn][oitr()]; oitr++;
#if BL_SPACEDIM>2
	const Mask& mb = *maskvals[level][gn][oitr()]; oitr++;
#endif
	const Mask& me = *maskvals[level][gn][oitr()]; oitr++;
	const Mask& mn = *maskvals[level][gn][oitr()]; oitr++;
#if BL_SPACEDIM>2
	const Mask& mt = *maskvals[level][gn][oitr()]; oitr++;
#endif
	FORT_DVFLUX(
	    xmfi().dataPtr(), 
	    ARLIM(xmfi().loVect()), ARLIM(xmfi().hiVect()),
	    &alpha, &beta,
	    amfi().dataPtr(), 
	    ARLIM(amfi().loVect()), ARLIM(amfi().hiVect()),
	    bXmfi().dataPtr(), 
	    ARLIM(bXmfi().loVect()), ARLIM(bXmfi().hiVect()),
	    bYmfi().dataPtr(), 
	    ARLIM(bYmfi().loVect()), ARLIM(bYmfi().hiVect()),
#if BL_SPACEDIM>2
	    bZmfi().dataPtr(), 
	    ARLIM(bZmfi().loVect()), ARLIM(bZmfi().hiVect()),
#endif
	    xfluxmfi().dataPtr(), 
	    ARLIM(xfluxmfi().loVect()), ARLIM(xfluxmfi().hiVect()),
	    yfluxmfi().dataPtr(), 
	    ARLIM(yfluxmfi().loVect()), ARLIM(yfluxmfi().hiVect()),
#if BL_SPACEDIM>2
	    zfluxmfi().dataPtr(), 
	    ARLIM(zfluxmfi().loVect()), ARLIM(zfluxmfi().hiVect()),
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
	    tdnfsi().dataPtr(),
	    ARLIM(tdnfsi().loVect()),ARLIM(tdnfsi().hiVect()),
	    tdefsi().dataPtr(),
	    ARLIM(tdefsi().loVect()),ARLIM(tdefsi().hiVect()),
	    tdwfsi().dataPtr(),
	    ARLIM(tdwfsi().loVect()),ARLIM(tdwfsi().hiVect()),
	    tdsfsi().dataPtr(),
	    ARLIM(tdsfsi().loVect()),ARLIM(tdsfsi().hiVect()),
#if BL_SPACEDIM>2
	    tdtfsi().dataPtr(),
	    ARLIM(tdtfsi().loVect()),ARLIM(tdtfsi().hiVect()),
	    tdbfsi().dataPtr(),
	    ARLIM(tdbfsi().loVect()),ARLIM(tdbfsi().hiVect()),
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
    const MultiFab& a  = aCoefficients(level);
    const MultiFab& bX = bCoefficients(0,level);
    const MultiFab& bY = bCoefficients(1,level);
#if BL_SPACEDIM>2
    const MultiFab& bZ = bCoefficients(2,level);
#endif
    OrientationIter oitr;
    const FabSet& tdw = (*tangderiv[level])[oitr()]; oitr++;
    const FabSet& tds = (*tangderiv[level])[oitr()]; oitr++;
#if BL_SPACEDIM>2
    const FabSet& tdb = (*tangderiv[level])[oitr()]; oitr++;
#endif
    const FabSet& tde = (*tangderiv[level])[oitr()]; oitr++;
    const FabSet& tdn = (*tangderiv[level])[oitr()]; oitr++;
#if BL_SPACEDIM>2
    const FabSet& tdt = (*tangderiv[level])[oitr()]; oitr++;
#endif
    const int nc = y.nComp();
    //
    // HACK: Cast away const for compatibility with BoxLib constructors.
    //
    for (MultiFabIterator xmfi((MultiFab&)x); xmfi.isValid(); ++xmfi)
    {
        DependentMultiFabIterator ymfi(xmfi,  y);
        DependentMultiFabIterator amfi(xmfi,  a);
        DependentMultiFabIterator bXmfi(xmfi, bX);
        DependentMultiFabIterator bYmfi(xmfi, bY);
#if BL_SPACEDIM > 2
        DependentMultiFabIterator bZmfi(xmfi, bZ);
#endif
	
        DependentFabSetIterator tdwfsi(xmfi, tdw);
        DependentFabSetIterator tdsfsi(xmfi, tds);
#if BL_SPACEDIM > 2
        DependentFabSetIterator tdbfsi(xmfi, tdb);
#endif
        DependentFabSetIterator tdefsi(xmfi, tde);
        DependentFabSetIterator tdnfsi(xmfi, tdn);
#if BL_SPACEDIM > 2
        DependentFabSetIterator tdtfsi(xmfi, tdt);
#endif

        oitr.rewind();
        const int gn = xmfi.index();
	const Mask& mw = *maskvals[level][gn][oitr()]; oitr++;
	const Mask& ms = *maskvals[level][gn][oitr()]; oitr++;
#if BL_SPACEDIM>2
	const Mask& mb = *maskvals[level][gn][oitr()]; oitr++;
#endif
	const Mask& me = *maskvals[level][gn][oitr()]; oitr++;
	const Mask& mn = *maskvals[level][gn][oitr()]; oitr++;
#if BL_SPACEDIM>2
	const Mask& mt = *maskvals[level][gn][oitr()]; oitr++;
#endif
	FORT_DVAPPLY(
	    xmfi().dataPtr(), 
            ARLIM(xmfi().loVect()), ARLIM(xmfi().hiVect()),
	    &alpha, &beta,
	    amfi().dataPtr(), 
            ARLIM(amfi().loVect()), ARLIM(amfi().hiVect()),
	    bXmfi().dataPtr(), 
            ARLIM(bXmfi().loVect()), ARLIM(bXmfi().hiVect()),
	    bYmfi().dataPtr(), 
            ARLIM(bYmfi().loVect()), ARLIM(bYmfi().hiVect()),
#if BL_SPACEDIM>2
	    bZmfi().dataPtr(), 
            ARLIM(bZmfi().loVect()), ARLIM(bZmfi().hiVect()),
#endif
	    ymfi().dataPtr(), 
            ARLIM(ymfi().loVect()), ARLIM(ymfi().hiVect()),
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
	    tdnfsi().dataPtr(),
	    ARLIM(tdnfsi().loVect()),ARLIM(tdnfsi().hiVect()),
	    tdefsi().dataPtr(),
	    ARLIM(tdefsi().loVect()),ARLIM(tdefsi().hiVect()),
	    tdwfsi().dataPtr(),
	    ARLIM(tdwfsi().loVect()),ARLIM(tdwfsi().hiVect()),
	    tdsfsi().dataPtr(),
	    ARLIM(tdsfsi().loVect()),ARLIM(tdsfsi().hiVect()),
#if BL_SPACEDIM>2
	    tdtfsi().dataPtr(),
	    ARLIM(tdtfsi().loVect()),ARLIM(tdtfsi().hiVect()),
	    tdbfsi().dataPtr(),
	    ARLIM(tdbfsi().loVect()),ARLIM(tdbfsi().hiVect()),
#endif
            xmfi.validbox().loVect(), xmfi.validbox().hiVect(),
	    h[level]
	    );
    }
}
