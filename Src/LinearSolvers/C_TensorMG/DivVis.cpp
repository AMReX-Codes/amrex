// $Id: DivVis.cpp,v 1.1 1998-03-26 19:15:32 marc Exp $

#include <DivVis.H>
#include <DivVis_F.H>

Real DivVis::a_def=0.0;
Real DivVis::b_def=1.0;
Real DivVis::alpha_def=1.0;
Real DivVis::beta_def=1.0;


int
DivVis::numberComponents()
{
  return BL_SPACEDIM;
}

int
DivVis::numberPhases()
{
#if BL_SPACEDIM==2
  return 4;
#else
  return 8;
#endif
}

DivVis::DivVis(const BoxArray &_ba, const MCBndryData &_bd,
                             Real _h)
    : MCLinOp(_ba, _bd, _h), alpha(alpha_def), beta(beta_def)
{
    initCoefficients(_ba);
    numcomp = numberComponents(); //wyc
    numphase = numberPhases(); // wyc
    undrrelxr.resize(1);
    undrrelxr[0] = new BndryRegister(*gbox[0], 1, 0, 0, numcomp);
    tangderiv.resize(1);
#if BL_SPACEDIM==2
    tangderiv[0] = new BndryRegister(*gbox[0], 0, 1, 0, numcomp);
#elif BL_SPACEDIM==3
    tangderiv[0] = new BndryRegister(*gbox[0], 0, 1, 0, numcomp*(1+3));
#else
# error
#endif
}

DivVis::DivVis(const BoxArray &_ba, const MCBndryData &_bd,
                             const Real* _h)
    : MCLinOp(_ba, _bd, _h), alpha(alpha_def), beta(beta_def)
{
    initCoefficients(_ba);
    numcomp = numberComponents(); //wyc
    numphase = numberPhases(); // wyc
    undrrelxr.resize(1);
    undrrelxr[0] = new BndryRegister(*gbox[0], 1, 0, 0, numcomp);
    tangderiv.resize(1);
#if BL_SPACEDIM==2
    tangderiv[0] = new BndryRegister(*gbox[0], 0, 1, 0, numcomp);
#elif BL_SPACEDIM==3
    tangderiv[0] = new BndryRegister(*gbox[0], 0, 1, 0, numcomp*(1+3));
#else
# error
#endif
}

DivVis::~DivVis()
{
    clearToLevel(-1);
}

void
DivVis::setScalars(Real _alpha, Real _beta)
{
    alpha = _alpha;
    beta  = _beta;
}

void
DivVis::clearToLevel(int level)
{
    assert(level>=-1);
    for (int i=level+1; i<numLevels(); ++i) {
	delete acoefs[i];
	for(int j=0; j<BL_SPACEDIM; ++j) {
	    delete bcoefs[i][j];
	}
    }
}

void
DivVis::prepareForLevel(int level)
{
    MCLinOp::prepareForLevel(level);
    if (level == 0 ) return;
    prepareForLevel(level-1);

      // If coefficients were marked invalid, or if not yet made, make new ones
      // (Note: makeCoefficients is a MCLinOp routine, and it allocates AND
      // fills coefficients.  A more efficient implementation would allocate
      // and fill in separate steps--we could then use the a_valid bool
      // along with the length of a_valid to separately determine whether to
      // fill or allocate the coefficient MultiFabs
    if (level >= a_valid.length() || a_valid[level] == false )  {
	if (acoefs.length() < level+1) {
	    acoefs.resize(level+1);
	    acoefs[level] = new MultiFab;
	} else {
	    delete acoefs[level];
	    acoefs[level] = new MultiFab;
	}
	makeCoefficients(*acoefs[level], *acoefs[level-1], level);
	a_valid.resize(level+1);
	a_valid[level] = true;
    }
    
    if (level >= b_valid.length() || b_valid[level] == false )  {
	if (bcoefs.length() < level+1) {
	    bcoefs.resize(level+1);
	    for(int i = 0; i < BL_SPACEDIM; ++i)
		bcoefs[level][i] = new MultiFab;
	} else {
	    for(int i = 0; i < BL_SPACEDIM; ++i) {
		delete bcoefs[level][i];
		bcoefs[level][i] = new MultiFab;
	    }
	}
	for(int i = 0; i < BL_SPACEDIM; ++i) {
	    makeCoefficients(*bcoefs[level][i], *bcoefs[level-1][i], level);
	}
	b_valid.resize(level+1);
	b_valid[level] = true;
    }
}

void
DivVis::initCoefficients(const BoxArray &_ba)
{
    int nGrow=0;
    acoefs.resize(1);
    bcoefs.resize(1);
    acoefs[0] = new MultiFab(_ba, 2, nGrow, Fab_allocate);
    acoefs[0]->setVal(a_def);
    a_valid.resize(1);
    a_valid[0] = true;

    for(int i = 0; i < BL_SPACEDIM; ++i) {
	BoxArray edge_boxes(_ba);
	edge_boxes.surroundingNodes(i);
	bcoefs[0][i] = new MultiFab(edge_boxes, 1, nGrow, Fab_allocate);
	bcoefs[0][i]->setVal(b_def);
    }
    b_valid.resize(1);
    b_valid[0] = true;
}

void
DivVis::invalidate_a_to_level(int lev)
{
    lev = (lev >= 0 ? lev : 0);
    for (int i=lev; i<numLevels(); i++) a_valid[i]=false;
}

void
DivVis::invalidate_b_to_level(int lev)
{
    lev = (lev >= 0 ? lev : 0);
    for (int i=lev; i<numLevels(); i++) b_valid[i]=false;
}


void
DivVis::aCoefficients(const MultiFab &_a)
{
    assert( _a.ok() );
    assert( _a.boxArray() == (acoefs[0])->boxArray() );
    invalidate_a_to_level(0);
    int ngrd = _a.length();
    for (int k = 0; k < ngrd; k++) {
      (*acoefs[0])[k].copy(_a[k],0,0,2);
    }
    return;
}

void
DivVis::bCoefficients(const MultiFab &_b, int dir)
{
    assert( _b.ok() );
    assert( _b.boxArray() == (bcoefs[0][dir])->boxArray() );
    invalidate_b_to_level(0);
    int ngrd = _b.length();
    for (int k = 0; k < ngrd; k++) {
      (*bcoefs[0][dir])[k].copy(_b[k],0,0,1);
    }
    return;
}

const MultiFab&
DivVis::aCoefficients(int level)
{
    prepareForLevel(level);
    return *acoefs[level];
}

const MultiFab&
DivVis::bCoefficients(int dir,int level)
{
    prepareForLevel(level);
    return *bcoefs[level][dir];
}

// must be defined for MultiGrid/CGSolver to work
void
DivVis::Fsmooth(MultiFab &solnL, const MultiFab &rhsL,
		       int level, int phaseflag)
{
    const BoxArray &bxa = *gbox[level];
    OrientationIter oitr;
    const FabSet &fw = (*undrrelxr[level])[oitr()]; 
    const FabSet &tdw = (*tangderiv[level])[oitr()];
    oitr++;
    const FabSet &fs = (*undrrelxr[level])[oitr()]; 
    const FabSet &tds = (*tangderiv[level])[oitr()];
    oitr++;
#if BL_SPACEDIM>2
    const FabSet &fb = (*undrrelxr[level])[oitr()]; 
    const FabSet &tdb = (*tangderiv[level])[oitr()];
    oitr++;
#endif
    const FabSet &fe = (*undrrelxr[level])[oitr()]; 
    const FabSet &tde = (*tangderiv[level])[oitr()];
    oitr++;
    const FabSet &fn = (*undrrelxr[level])[oitr()]; 
    const FabSet &tdn = (*tangderiv[level])[oitr()];
    oitr++;
#if BL_SPACEDIM>2
    const FabSet &ft = (*undrrelxr[level])[oitr()]; 
    const FabSet &tdt = (*tangderiv[level])[oitr()];
    oitr++;
#endif
    const MultiFab  &a = aCoefficients(level);
    const MultiFab  &bX = bCoefficients(0,level);
    const MultiFab  &bY = bCoefficients(1,level);
#if BL_SPACEDIM>2
    const MultiFab  &bZ = bCoefficients(2,level);
#endif
    int nc = solnL.nComp();
    for(int gn = 0; gn < solnL.length(); ++gn) {
	oitr.rewind();
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

	FArrayBox &s = solnL[gn];
	Real* sptr = s.dataPtr();
	const Real* rhsptr = rhsL[gn].dataPtr(); 
	const Real* aptr = a[gn].dataPtr();  
	const Real* bXptr = bX[gn].dataPtr(); 
	const Real* bYptr = bY[gn].dataPtr(); 
#if BL_SPACEDIM>2
	const Real* bZptr = bZ[gn].dataPtr(); 
#endif

	const Real* fnptr = fn[gn].dataPtr();
	const int* mnptr = mn.dataPtr(); 
	const Real* feptr = fe[gn].dataPtr();
	const int* meptr = me.dataPtr(); 
	const Real* fwptr = fw[gn].dataPtr();
	const int* mwptr = mw.dataPtr(); 
	const Real* fsptr = fs[gn].dataPtr();
	const int* msptr = ms.dataPtr(); 
#if BL_SPACEDIM>2
	const Real* ftptr = ft[gn].dataPtr();
	const int* mtptr = mt.dataPtr(); 
	const Real* fbptr = fb[gn].dataPtr();
	const int* mbptr = mb.dataPtr(); 
#endif

	const Real* tdnptr = tdn[gn].dataPtr();
	const FArrayBox& tefab = tde[gn];
	const Real* tdeptr = tefab.dataPtr();
	const FArrayBox& twfab = tdw[gn];
	const Real* tdwptr = twfab.dataPtr();
	const Real* tdsptr = tds[gn].dataPtr();
#if BL_SPACEDIM>2
	const Real* tdtptr = tdt[gn].dataPtr();
	const Real* tdbptr = tdb[gn].dataPtr();
#endif

	FORT_GSRB(
	    sptr, 
            ARLIM(s.loVect()),ARLIM(s.hiVect()),
	    rhsptr,
            ARLIM(rhsL[gn].loVect()), ARLIM(rhsL[gn].hiVect()),
	    &alpha, &beta,
	    aptr,
            ARLIM(a[gn].loVect()),    ARLIM(a[gn].hiVect()),
	    bXptr,
            ARLIM(bX[gn].loVect()),   ARLIM(bX[gn].hiVect()),
	    bYptr,
            ARLIM(bY[gn].loVect()),   ARLIM(bY[gn].hiVect()),
#if BL_SPACEDIM>2
	    bZptr,
            ARLIM(bZ[gn].loVect()),   ARLIM(bZ[gn].hiVect()),
#endif
	    mnptr,
            ARLIM(mn.loVect()),   ARLIM(mn.hiVect()),
	    fnptr,
            ARLIM(fn[gn].loVect()),   ARLIM(fn[gn].hiVect()),
	    meptr,
            ARLIM(me.loVect()),   ARLIM(me.hiVect()),
	    feptr,
            ARLIM(fe[gn].loVect()),   ARLIM(fe[gn].hiVect()),
	    mwptr,
            ARLIM(mw.loVect()),   ARLIM(mw.hiVect()),
	    fwptr,
            ARLIM(fw[gn].loVect()),   ARLIM(fw[gn].hiVect()),
	    msptr,
            ARLIM(ms.loVect()),   ARLIM(ms.hiVect()),
	    fsptr,
            ARLIM(fs[gn].loVect()),   ARLIM(fs[gn].hiVect()),
#if BL_SPACEDIM>2
	    mtptr,
            ARLIM(mt.loVect()),   ARLIM(mt.hiVect()),
	    ftptr,
            ARLIM(ft[gn].loVect()),   ARLIM(ft[gn].hiVect()),
	    mbptr,
            ARLIM(mb.loVect()),   ARLIM(mb.hiVect()),
	    fbptr,
            ARLIM(fb[gn].loVect()),   ARLIM(fb[gn].hiVect()),
#endif
	    tdnptr,
	    ARLIM(tdn[gn].loVect()),ARLIM(tdn[gn].hiVect()),
	    tdeptr,
	    ARLIM(tde[gn].loVect()),ARLIM(tde[gn].hiVect()),
	    tdwptr,
	    ARLIM(tdw[gn].loVect()),ARLIM(tdw[gn].hiVect()),
	    tdsptr,
	    ARLIM(tds[gn].loVect()),ARLIM(tds[gn].hiVect()),
#if BL_SPACEDIM>2
	    tdtptr,
	    ARLIM(tdt[gn].loVect()),ARLIM(tdt[gn].hiVect()),
	    tdbptr,
	    ARLIM(tdb[gn].loVect()),ARLIM(tdb[gn].hiVect()),
#endif
	    bxa[gn].loVect(), bxa[gn].hiVect(), h[level], nc, phaseflag);

    }
}

void 
DivVis::compFlux(MultiFab &xflux, 
		 MultiFab &yflux, 
#if BL_SPACEDIM>2
		 MultiFab &zflux, 
#endif
		 MultiFab& x)
{
  int level = 0;
  MCBC_Mode bc_mode = MCInhomogeneous_BC;
  applyBC(x,level,bc_mode);

  const BoxArray &bxa = *gbox[level];
  const MultiFab  &a = aCoefficients(level);
  const MultiFab  &bX = bCoefficients(0,level);
  const MultiFab  &bY = bCoefficients(1,level);
#if BL_SPACEDIM>2
  const MultiFab  &bZ = bCoefficients(2,level);
#endif
  OrientationIter oitr;
  const FabSet &tdw = (*tangderiv[level])[oitr()]; oitr++;
  const FabSet &tds = (*tangderiv[level])[oitr()]; oitr++;
#if BL_SPACEDIM>2
  const FabSet &tdb = (*tangderiv[level])[oitr()]; oitr++;
#endif
  const FabSet &tde = (*tangderiv[level])[oitr()]; oitr++;
  const FabSet &tdn = (*tangderiv[level])[oitr()]; oitr++;
#if BL_SPACEDIM>2
  const FabSet &tdt = (*tangderiv[level])[oitr()]; oitr++;
#endif
  int nc = x.nComp();
  assert( nc == BL_SPACEDIM );
  assert( nc == xflux.nComp() );
  assert( nc == yflux.nComp() );

  for(int gn = 0; gn < x.length(); ++gn) {
    oitr.rewind();
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
		x[gn].dataPtr(), 
		ARLIM(x[gn].loVect()), ARLIM(x[gn].hiVect()),
		&alpha, &beta,
		a[gn].dataPtr(), 
		ARLIM(a[gn].loVect()), ARLIM(a[gn].hiVect()),
		bX[gn].dataPtr(), 
		ARLIM(bX[gn].loVect()), ARLIM(bX[gn].hiVect()),
		bY[gn].dataPtr(), 
		ARLIM(bY[gn].loVect()), ARLIM(bY[gn].hiVect()),
#if BL_SPACEDIM>2
		bZ[gn].dataPtr(), 
		ARLIM(bZ[gn].loVect()), ARLIM(bZ[gn].hiVect()),
#endif
		xflux[gn].dataPtr(), 
		ARLIM(xflux[gn].loVect()), ARLIM(xflux[gn].hiVect()),
		yflux[gn].dataPtr(), 
		ARLIM(yflux[gn].loVect()), ARLIM(yflux[gn].hiVect()),
#if BL_SPACEDIM>2
		yflux[gn].dataPtr(), 
		ARLIM(yflux[gn].loVect()), ARLIM(yflux[gn].hiVect()),
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
		tdn[gn].dataPtr(),
		ARLIM(tdn[gn].loVect()),ARLIM(tdn[gn].hiVect()),
		tde[gn].dataPtr(),
		ARLIM(tde[gn].loVect()),ARLIM(tde[gn].hiVect()),
		tdw[gn].dataPtr(),
		ARLIM(tdw[gn].loVect()),ARLIM(tdw[gn].hiVect()),
		tds[gn].dataPtr(),
		ARLIM(tds[gn].loVect()),ARLIM(tds[gn].hiVect()),
#if BL_SPACEDIM>2
		tdt[gn].dataPtr(),
		ARLIM(tdt[gn].loVect()),ARLIM(tdt[gn].hiVect()),
		tdb[gn].dataPtr(),
		ARLIM(tdb[gn].loVect()),ARLIM(tdb[gn].hiVect()),
#endif
		bxa[gn].loVect(), bxa[gn].hiVect(), 
		h[level]
		);
  }
}



void
DivVis::Fapply(MultiFab& y, const MultiFab &x, int level)
{
    const BoxArray &bxa = *gbox[level];
    const MultiFab  &a = aCoefficients(level);
    const MultiFab  &bX = bCoefficients(0,level);
    const MultiFab  &bY = bCoefficients(1,level);
#if BL_SPACEDIM>2
    const MultiFab  &bZ = bCoefficients(2,level);
#endif
    OrientationIter oitr;
    const FabSet &tdw = (*tangderiv[level])[oitr()]; oitr++;
    const FabSet &tds = (*tangderiv[level])[oitr()]; oitr++;
#if BL_SPACEDIM>2
    const FabSet &tdb = (*tangderiv[level])[oitr()]; oitr++;
#endif
    const FabSet &tde = (*tangderiv[level])[oitr()]; oitr++;
    const FabSet &tdn = (*tangderiv[level])[oitr()]; oitr++;
#if BL_SPACEDIM>2
    const FabSet &tdt = (*tangderiv[level])[oitr()]; oitr++;
#endif
    int nc = y.nComp();
    for(int gn = 0; gn < y.length(); ++gn) {
      oitr.rewind();
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
	    x[gn].dataPtr(), 
            ARLIM(x[gn].loVect()), ARLIM(x[gn].hiVect()),
	    &alpha, &beta,
	    a[gn].dataPtr(), 
            ARLIM(a[gn].loVect()), ARLIM(a[gn].hiVect()),
	    bX[gn].dataPtr(), 
            ARLIM(bX[gn].loVect()), ARLIM(bX[gn].hiVect()),
	    bY[gn].dataPtr(), 
            ARLIM(bY[gn].loVect()), ARLIM(bY[gn].hiVect()),
#if BL_SPACEDIM>2
	    bZ[gn].dataPtr(), 
            ARLIM(bZ[gn].loVect()), ARLIM(bZ[gn].hiVect()),
#endif
	    y[gn].dataPtr(), 
            ARLIM(y[gn].loVect()), ARLIM(y[gn].hiVect()),
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
	    tdn[gn].dataPtr(),
	    ARLIM(tdn[gn].loVect()),ARLIM(tdn[gn].hiVect()),
	    tde[gn].dataPtr(),
	    ARLIM(tde[gn].loVect()),ARLIM(tde[gn].hiVect()),
	    tdw[gn].dataPtr(),
	    ARLIM(tdw[gn].loVect()),ARLIM(tdw[gn].hiVect()),
	    tds[gn].dataPtr(),
	    ARLIM(tds[gn].loVect()),ARLIM(tds[gn].hiVect()),
#if BL_SPACEDIM>2
	    tdt[gn].dataPtr(),
	    ARLIM(tdt[gn].loVect()),ARLIM(tdt[gn].hiVect()),
	    tdb[gn].dataPtr(),
	    ARLIM(tdb[gn].loVect()),ARLIM(tdb[gn].hiVect()),
#endif
	    bxa[gn].loVect(), bxa[gn].hiVect(), 
	    h[level]
	    );
    }
}
