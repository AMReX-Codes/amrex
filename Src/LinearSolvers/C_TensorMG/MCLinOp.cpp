// $Id: MCLinOp.cpp,v 1.1 1998-03-26 19:15:36 marc Exp $

// differences from LinOp: den has nc components, bct has nc components

#include <stdlib.h>
#include <ParmParse.H>

#include "LO_BCTYPES.H"
#include "DivVis_F.H"
#include "MCLinOp.H"

bool MCLinOp::initialized = false;
int MCLinOp::def_harmavg = 0;
int MCLinOp::def_verbose = 0;
int MCLinOp::def_maxorder = 2;

  // MCLinOp::applyBC fills MCLinOp_grow ghost cells with data expected in MCLinOp::apply
  //  therefore, the incoming MultiFab to MCLinOp::applyBC better have this many ghost
  //  allocated
const int MCLinOp_grow = 1;

  // MCLinOp::computePeriodicIntersections precomputes some intersection information
  // assuming that Periodic_BC_grow cells are important enough to be included in
  // the fix-up after the MultiFab::filBndry() call in MCLinOp::applyBC.  Since
  // MCLinOp::applyBC is presently the only member requiring ghostcells, Periodic_BC_grow
  // is set to its maximum requirement.
const int Periodic_BC_grow = 1;

void
MCLinOp::initialize()
{
    ParmParse pp("MCLp");
    pp.query("harmavg", def_harmavg);
    pp.query("verbose", def_verbose);
    pp.query("v", def_verbose);
    if (def_verbose) {
	if(ParallelDescriptor::IOProcessor()) {
	    cout << "def_harmavg = " << def_harmavg << '\n';
	}
    }
    initialized = true;
}

MCLinOp::MCLinOp(const BoxArray &ba, const MCBndryData &_bgb, const Real _h)
    : bgb(_bgb)
{
    if(!initialized)
	initialize();
    harmavg = def_harmavg;
    verbose = def_verbose;
    gbox.resize(1);
    gbox[0] = new BoxArray(ba);
    geomarray.resize(1);
    geomarray[0] = bgb.getGeom();
    pirarray.resize(1);
    pirarray[0].resize(100);  // make some extra room
    pirarray[0].resize(0);  // make some extra room
    computePeriodicIntersections( geomarray[0],*(gbox[0]),pirarray[0] );
    h.resize(1);
    maxorder = def_maxorder;
    int i;
    for(i = 0; i < BL_SPACEDIM; ++i) {
	h[0][i] = _h;
    }
    maskvals.resize(1);
    maskvals[0].resize(bgb.boxes().length());
    for(i=0; i < bgb.boxes().length(); ++i) {
	OrientationIter oitr;
	maskvals[0][i].resize(2*BL_SPACEDIM, (Mask*)0);
	int m = 0;
	while ( oitr ) {
	    const PArray<Mask> &pam = bgb.bndryMasks(oitr());
	    maskvals[0][i][m] = new Mask(pam[i].box(), 1);
	    maskvals[0][i][m++]->copy(pam[i]);
	    oitr++;
	}
    }
}

MCLinOp::MCLinOp(const BoxArray &ba, const MCBndryData &_bgb, const Real* _h)
    : bgb(_bgb)
{
    if(!initialized)
	initialize();
    harmavg = def_harmavg;
    verbose = def_verbose;
    gbox.resize(1);
    gbox[0] = new BoxArray(ba);
    geomarray.resize(1);
    geomarray[0] = bgb.getGeom();
    pirarray.resize(1);
    pirarray[0].resize(100);  // make some extra room
    pirarray[0].resize(0);  // make some extra room
    computePeriodicIntersections( geomarray[0],*(gbox[0]),pirarray[0] );
    h.resize(1);
    maxorder = def_maxorder;
    int i;
    for(i = 0; i < BL_SPACEDIM; ++i) {
	h[0][i] = _h[i];
    }
    maskvals.resize(1);
    maskvals[0].resize(bgb.boxes().length());
    for(i=0; i < bgb.boxes().length(); ++i) {
	OrientationIter oitr;
	maskvals[0][i].resize(2*BL_SPACEDIM, (Mask*)0);
	int m = 0;
	while ( oitr ) {
	    const PArray<Mask> &pam = bgb.bndryMasks(oitr());
	    maskvals[0][i][m] = new Mask(pam[i].box(), 1);
	    maskvals[0][i][m++]->copy(pam[i]);
	    oitr++;
	}
    }
}

MCLinOp::~MCLinOp()
{
    int i;
    for(i=0; i < maskvals.length(); ++i) {
	for(int j=0; j < maskvals[i].length(); ++j) {
	    for(int k = 0; k < maskvals[i][j].length(); ++k) {
		delete maskvals[i][j][k];
	    }
	}
    }
}

MCLinOp::MCLinOp(const MCLinOp& _lp, int level)
    : bgb(_lp.bgb)
{
    harmavg = _lp.harmavg;
    verbose = _lp.verbose;
    gbox.resize(1);
    gbox[0] = new BoxArray(_lp.boxArray(level));
    geomarray.resize(1);
    geomarray[0] = bgb.getGeom();
    pirarray.resize(1);
    pirarray[0].resize(100);  // make some extra room
    pirarray[0].resize(0);  // make some extra room
    h.resize(1);
    h[0] = _lp.h[level];	// level should be prepared here.
    undrrelxr.resize(1);
    undrrelxr[0] = _lp.undrrelxr[level];
    tangderiv.resize(1);
    tangderiv[0] = _lp.tangderiv[level];
}

void
MCLinOp::apply(MultiFab& out, MultiFab& in, int level, MCBC_Mode bc_mode)
{
    applyBC(in,level,bc_mode);
    Fapply(out,in,level);
}

void
MCLinOp::applyBC(MultiFab &inout, int level, MCBC_Mode bc_mode)
{
      // The inout MultiFab must have at least MCLinOp_grow ghost cells for applyBC
    assert( inout.nGrow() >= MCLinOp_grow);
    
      // The inout MultiFab must have at least Periodic_BC_grow cells for the
      // algorithms taking care of periodic boundary conditions
    assert( inout.nGrow() >= MCLinOp_grow);
    
      // No coarsened boundary values, cannot apply inhomog at lev>0
    assert( !(level>0 && bc_mode == MCInhomogeneous_BC) );
    
    int flagden = 1;	// fill in the bndry data and undrrelxr
    int flagbc  = 1;	// with values
    if (bc_mode == MCHomogeneous_BC) flagbc = 0; // nodata if homog
    int nc = inout.nComp();
    assert(nc == numcomp );

    inout.FillBoundary();
    prepareForLevel(level);

    // do periodic fixup
    const Geometry& geom = geomarray[level];
    if( geom.isAnyPeriodic() ){
      ParallelDescriptor::Abort("LinOp::applyBC not yet implemented for periodic in parallel");
      const Array<struct PeriodicIntersectionRec> & pir = pirarray[level];
      int len = pir.length();
      for( int i=0; i<len; i++){
	const struct PeriodicIntersectionRec &rec = pir[i];
	const Box & srcbox = rec.srcbox;
	const Box & destbox = rec.destbox;
	int srcno = rec.srcno;
	int destno = rec.destno;
	inout[destno].copy(inout[srcno],srcbox,0,destbox,0,nc);
      }
    }

    OrientationIter oitr;
    while( oitr ) {
	const Array<Array<BoundCond> > &b = bgb.bndryConds(oitr());
	const Array<Real> &r = bgb.bndryLocs(oitr());
	FabSet &f = (*undrrelxr[level])[oitr()];
	FabSet &td = (*tangderiv[level])[oitr()];
	int cdr(oitr());
	const FabSet& fs = bgb.bndryValues(oitr());
	int cdir = oitr().coordDir();
        for(MultiFabIterator inoutmfi(inout); inoutmfi.isValid(); ++inoutmfi) {
          DependentFabSetIterator ffsi(inoutmfi, f);
          DependentFabSetIterator tdfsi(inoutmfi, td);
          DependentFabSetIterator fsfsi(inoutmfi, fs);
	    int gn = inoutmfi.index();
	    assert(gbox[level]->get(inoutmfi.index()) == inoutmfi.validbox());
	    Real bcl(r[gn]);
	    const int *bct = (const int *)b[gn].dataPtr();
	    FArrayBox& fsfab = fsfsi();
	    const Real* bcvalptr = fsfab.dataPtr();
	    // way external derivs stored
	    const Real* exttdptr = fsfab.dataPtr(numcomp); 
	    const int* fslo = fsfab.loVect();
	    const int* fshi = fsfab.hiVect();
	    FArrayBox& inoutfab = inoutmfi();
	    FArrayBox& denfab = ffsi();
	    FArrayBox& tdfab = tdfsi();
#if BL_SPACEDIM==2
	    int perpdir;
	    if( cdir == 0 ) perpdir = 1;
	    else if( cdir == 1 ) perpdir = 0;
	    else {
	      cerr << "MCLinOp::applyBC: bad logic"<<endl;
	      exit(1);
	    }
	    const Mask& m = *maskvals[level][gn][oitr()];
	    const Mask& mphi = *maskvals[level][gn][Orientation(perpdir,
							Orientation::high)];
	    const Mask& mplo = *maskvals[level][gn][Orientation(perpdir,
							Orientation::low)];
	    FORT_APPLYBC(
		&flagden, &flagbc, &maxorder,
		inoutfab.dataPtr(), 
                ARLIM(inoutfab.loVect()), ARLIM(inoutfab.hiVect()),
		&cdr, bct, &bcl,
		bcvalptr, ARLIM(fslo), ARLIM(fshi),
		m.dataPtr(), ARLIM(m.loVect()), ARLIM(m.hiVect()),
		mphi.dataPtr(), ARLIM(mphi.loVect()), ARLIM(mphi.hiVect()),
		mplo.dataPtr(), ARLIM(mplo.loVect()), ARLIM(mplo.hiVect()),
		denfab.dataPtr(), 
		ARLIM(denfab.loVect()), ARLIM(denfab.hiVect()),
		exttdptr, ARLIM(fslo), ARLIM(fshi),
		tdfab.dataPtr(),ARLIM(tdfab.loVect()),ARLIM(tdfab.hiVect()),
		gbox[level]->get(gn).loVect(), gbox[level]->get(gn).hiVect(),
		&nc, h[level]);
#elif BL_SPACEDIM==3
	    const Mask& mn = *maskvals[level][gn][Orientation(1,
							 Orientation::high)];
	    const Mask& me = *maskvals[level][gn][Orientation(0,
							 Orientation::high)];
	    const Mask& mw = *maskvals[level][gn][Orientation(0,
							 Orientation::low)];
	    const Mask& ms = *maskvals[level][gn][Orientation(1,
							 Orientation::low)];
	    const Mask& mt = *maskvals[level][gn][Orientation(2,
							 Orientation::high)];
	    const Mask& mb = *maskvals[level][gn][Orientation(2,
							 Orientation::low)];
	    FORT_APPLYBC(
		&flagden, &flagbc, &maxorder,
		inoutfab.dataPtr(), 
                ARLIM(inoutfab.loVect()), ARLIM(inoutfab.hiVect()),
		&cdr, bct, &bcl,
		bcvalptr, ARLIM(fslo), ARLIM(fshi),
		mn.dataPtr(),ARLIM(mn.loVect()),ARLIM(mn.hiVect()),
		me.dataPtr(),ARLIM(me.loVect()),ARLIM(me.hiVect()),
		mw.dataPtr(),ARLIM(mw.loVect()),ARLIM(mw.hiVect()),
		ms.dataPtr(),ARLIM(ms.loVect()),ARLIM(ms.hiVect()),
		mt.dataPtr(),ARLIM(mt.loVect()),ARLIM(mt.hiVect()),
		mb.dataPtr(),ARLIM(mb.loVect()),ARLIM(mb.hiVect()),
		denfab.dataPtr(), 
		ARLIM(denfab.loVect()), ARLIM(denfab.hiVect()),
		exttdptr, ARLIM(fslo), ARLIM(fshi),
		tdfab.dataPtr(),ARLIM(tdfab.loVect()),ARLIM(tdfab.hiVect()),
		gbox[level]->get(gn).loVect(), gbox[level]->get(gn).hiVect(),
		&nc, h[level]);
#endif
	}
	oitr++;
    }
}
    
void
MCLinOp::residual(MultiFab &residL, const MultiFab &rhsL, MultiFab &solnL,
		int level, MCBC_Mode bc_mode)
{
    apply(residL, solnL, level, bc_mode);
    for(MultiFabIterator solnLmfi(solnL); solnLmfi.isValid(); ++solnLmfi) {
      DependentMultiFabIterator residLmfi(solnLmfi, residL);
      DependentMultiFabIterator rhsLmfi(solnLmfi, rhsL);
	int nc = residL.nComp();
	FORT_RESIDL(
	    residLmfi().dataPtr(), 
            ARLIM(residLmfi().loVect()), ARLIM(residLmfi().hiVect()),
	    rhsLmfi().dataPtr(), 
            ARLIM(rhsLmfi().loVect()), ARLIM(rhsLmfi().hiVect()),
	    residLmfi().dataPtr(), 
            ARLIM(residLmfi().loVect()), ARLIM(residLmfi().hiVect()),
	    solnLmfi.validbox().loVect(), solnLmfi.validbox().hiVect(), &nc
	    );
    }
}

void
MCLinOp::smooth(MultiFab &solnL, const MultiFab &rhsL,
	      int level, MCBC_Mode bc_mode)
{
    for (int phaseflag = 0; phaseflag < numphase; phaseflag++) {
      applyBC(solnL, level, bc_mode);
      Fsmooth(solnL, rhsL, level, phaseflag);
   }
}

Real
MCLinOp::norm(const MultiFab &in, int level) const
{
    Real norm = 0.0;
    for(ConstMultiFabIterator inmfi(in); inmfi.isValid(); ++inmfi) {
        int gn = inmfi.index();
        Real tnorm = inmfi().norm(gbox[level]->get(gn));
	norm += tnorm*tnorm;
    }
    return norm;
}

void
MCLinOp::clearToLevel(int level)
{
    int i;
    for(i = level+1; i < numLevels(); ++i) {
	delete undrrelxr[i].release();
	delete tangderiv[i].release();
	gbox[i] = 0;
    }
    h.resize(level+1);
    gbox.resize(level+1);
    undrrelxr.resize(level+1);
    tangderiv.resize(level+1);
}

void
MCLinOp::computePeriodicIntersections(const Geometry&geo, const BoxArray& boxarr,
				    Array<struct PeriodicIntersectionRec>&pir)
{
  if( ! geo.isAnyPeriodic() ) return;
  pir.resize(0);

  const Box& domain = geo.Domain();
  int len = boxarr.length();

  for( int i=0; i<len; i++ ){
    Box dest = boxarr[i];
    dest.grow(Periodic_BC_grow);
    if( ! domain.contains(dest) ){
      for( int j=0; j<len; j++ ){
	Box src = boxarr[j];
	Array<IntVect> pshifts;
	geo.periodicShift( dest, src, pshifts );
	for( int iiv=0; iiv<pshifts.length(); iiv++ ){
	  IntVect iv = pshifts[iiv];
	  Box shbox( src );
	  D_TERM( shbox.shift(0,iv[0]);,
		  shbox.shift(1,iv[1]);,
		  shbox.shift(2,iv[2]); )
	  Box intbox = dest & shbox;
	  assert( intbox.ok() );
	  // ok, we got an intersection
	  pir.resize(pir.length()+1);
	  struct PeriodicIntersectionRec & rec = pir[pir.length()-1];
	  rec.srcno = j;
	  rec.destno = i;
	  rec.destbox = intbox;
	  D_TERM( intbox.shift(0,-iv[0]);,
		  intbox.shift(1,-iv[1]);,
		  intbox.shift(2,-iv[2]); )
	  rec.srcbox = intbox;
	}
      }
    }
  }
}

void
MCLinOp::prepareForLevel(int level)
{
    if(level == 0)
	return;
    MCLinOp::prepareForLevel(level-1);
    if(h.length() > level)
	return;
    h.resize(level+1);
    int i;
    for(i = 0; i < BL_SPACEDIM; ++i) {
	h[level][i] = h[level-1][i]*2.0;
    }
    geomarray.resize(level+1);
    Box curdomain( geomarray[level-1].Domain() );
    curdomain.coarsen(2);
    geomarray[level].define( curdomain );
    gbox.resize(level+1);
    gbox[level] = new BoxArray(*(gbox[level-1]));
    gbox[level]->coarsen(2);

    pirarray.resize(level+1);
    pirarray[level].resize(100);  // make some extra room
    pirarray[level].resize(0);  // make some extra room
    computePeriodicIntersections(geomarray[level],*(gbox[level]),
				 pirarray[level]);

    undrrelxr.resize(level+1);
    undrrelxr[level] = new BndryRegister(*gbox[level], 1, 0, 0, numcomp);
    tangderiv.resize(level+1);
    // figure out how many components
    const FabSet& samplefs = (*tangderiv[level-1])[Orientation(0,
							    Orientation::low)];
    const FArrayBox & samplefab = samplefs[0];
    int ntdcomp = samplefab.nComp();

    tangderiv[level] = new BndryRegister(*gbox[level], 0, 1, 0, ntdcomp);
    maskvals.resize(level+1);
    maskvals[level].resize(gbox[level]->length());
    for(i=0; i < gbox[level]->length(); ++i) {
	maskvals[level][i].resize(2*BL_SPACEDIM, (Mask*)0);
    	int m = 0;
    	OrientationIter oitr;
    	while ( oitr ) {
	    BOX bx_k = adjCell((*gbox[level])[i], oitr(), 1);
	    // extend box in directions orthogonal to face normal
	    for( int dir=0; dir < BL_SPACEDIM; dir++ ){
	      if( dir==oitr() ) continue;
	      bx_k.grow(dir,1);
	    }
    	    maskvals[level][i][m] = new Mask(bx_k, 1);
    	    maskvals[level][i][m]->setVal(MCBndryData::not_covered);
	    for(int gn = 0; gn < gbox[level]->length(); ++gn) {
		BOX btmp = (*gbox[level])[gn];
		btmp &= bx_k;
		maskvals[level][i][m]->setVal(MCBndryData::covered, btmp,0);
	    }
	    // now take care of periodic wraparounds
	    Geometry& curgeom = geomarray[level];
	    if( curgeom.isAnyPeriodic() && !curdomain.contains(bx_k)  ){
	      Array<IntVect> pshifts(27);
	      curgeom.periodicShift( curdomain,bx_k,pshifts);
	      Mask &curmask = *(maskvals[level][i][m]);
	      for( int iiv=0; iiv<pshifts.length(); iiv++ ){
		IntVect iv = pshifts[iiv];
		curmask.shift(iv);

		for(int gn=0; gn<gbox[level]->length(); ++gn){
		  BOX btmp = (*gbox[level])[gn];
		  btmp &= curmask.box();
		  curmask.setVal(MCBndryData::covered, btmp,0);
		}
		
		curmask.shift(-iv);
	      }
	    }
	    oitr++;
	    m++;
	}
    }
}

void
MCLinOp::makeCoefficients(MultiFab& cs, const MultiFab &fn, int level)
{
    int nc = fn.nComp();
    const BoxArray& bxa = fn.boxArray();

      // Determine index type of incoming MultiFab
    IndexType itype = bxa[0].ixType();
    int cdir = -1;		// default to CELLTYPE
#if (BL_SPACEDIM == 2)
    if ( itype == IndexType(IndexType::NODE, IndexType::CELL) ) {
	cdir =  0;
    } else if ( itype == IndexType(IndexType::CELL, IndexType::NODE) ) {
	cdir =  1;
    }
#endif
#if (BL_SPACEDIM == 3)
    if (        itype == IndexType(IndexType::NODE,
				   IndexType::CELL,
				   IndexType::CELL) ) {
	cdir =  0;
    } else if ( itype == IndexType(IndexType::CELL,
				   IndexType::NODE,
				   IndexType::CELL) ) {
	cdir =  1;
    } else if ( itype == IndexType(IndexType::CELL,
				   IndexType::CELL,
				   IndexType::NODE) ) {
	cdir =  2;
    }
#endif
    
    BoxArray d(*gbox[level]);
    if(cdir >= 0) {
	d.surroundingNodes(cdir);
    }

    int nGrow=0;
    cs.define(d, nc, nGrow, Fab_allocate);
    cs.setVal(0.0);

      // some abbreviations...
    BoxArray &mgbl = *(gbox[level]);
    for(MultiFabIterator csmfi(cs); csmfi.isValid(); ++csmfi) {
      DependentMultiFabIterator fnmfi(csmfi, fn);
	switch(cdir) {
	case -1:
	    FORT_AVERAGECC(
		csmfi().dataPtr(), ARLIM(csmfi().loVect()), ARLIM(csmfi().hiVect()),
		fnmfi().dataPtr(), ARLIM(fnmfi().loVect()), ARLIM(fnmfi().hiVect()),
		mgbl[csmfi.index()].loVect(), mgbl[csmfi.index()].hiVect(), &nc
		);
	    break;
	case 0:
	case 1:
	case 2:
	    if ( harmavg ) {
		FORT_HARMONIC_AVERAGEEC(
		    csmfi().dataPtr(), 
                    ARLIM(csmfi().loVect()), ARLIM(csmfi().hiVect()),
		    fnmfi().dataPtr(), 
                    ARLIM(fnmfi().loVect()), ARLIM(fnmfi().hiVect()),
		    mgbl[csmfi.index()].loVect(), mgbl[csmfi.index()].hiVect(), &nc, &cdir
		    );
	    } else {
		FORT_AVERAGEEC(
		    csmfi().dataPtr(), 
                    ARLIM(csmfi().loVect()), ARLIM(csmfi().hiVect()),
		    fnmfi().dataPtr(), 
                    ARLIM(fnmfi().loVect()), ARLIM(fnmfi().hiVect()),
		    mgbl[csmfi.index()].loVect(), mgbl[csmfi.index()].hiVect(), &nc, &cdir
		    );
	    }
	    break;
	default:
	    BoxLib::Error("LinOp:: bad coefficient coarsening direction!");
	}
    }
}
