//
// $Id: LinOp.cpp,v 1.1 1998-03-24 07:05:35 almgren Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#include <ParmParse.H>
#include <ParallelDescriptor.H>

#include <LO_BCTYPES.H>
#include <LO_F.H>
#include <LinOp.H>

bool LinOp::initialized = false;
int LinOp::def_harmavg = 0;
int LinOp::def_verbose = 0;
int LinOp::def_maxorder = 2;

#ifndef NDEBUG
  // LinOp::applyBC fills LinOp_grow ghost cells with data expected in LinOp::apply
  //  therefore, the incoming MultiFab to LinOp::applyBC better have this many ghost
  //  allocated
const int LinOp_grow = 1;
#endif

  // LinOp::computePeriodicIntersections precomputes some intersection information
  // assuming that Periodic_BC_grow cells are important enough to be included in
  // the fix-up after the MultiFab::FillBoundary() call in LinOp::applyBC.  Since
  // LinOp::applyBC is presently the only member requiring ghostcells, Periodic_BC_grow
  // is set to its maximum requirement.
const int Periodic_BC_grow = 1;

void
LinOp::initialize()
{
    ParmParse pp("Lp");
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

LinOp::LinOp(const BoxArray &ba, const BndryData &_bgb, const Real _h)
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
    undrrelxr.resize(1);
    undrrelxr[0] = new BndryRegister(*gbox[0], 1, 0, 0, 1);
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

LinOp::LinOp(const BoxArray &ba, const BndryData &_bgb, const Real* _h)
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
    undrrelxr.resize(1);
    undrrelxr[0] = new BndryRegister(*gbox[0], 1, 0, 0, 1);
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

LinOp::~LinOp()
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

LinOp::LinOp(const LinOp& _lp, int level)
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
    h[0] = _lp.h[level];        // level should be prepared here.
    undrrelxr.resize(1);
    undrrelxr[0] = _lp.undrrelxr[level];
}

void
LinOp::apply(MultiFab& out, MultiFab& in, int level, LinOp::BC_Mode bc_mode)
{
    applyBC(in,level,bc_mode);
    Fapply(out,in,level);
}

void
LinOp::applyBC(MultiFab &inout, int level, LinOp::BC_Mode bc_mode)
{
      // The inout MultiFab must have at least LinOp_grow ghost cells for applyBC
    assert( inout.nGrow() >= LinOp_grow);
    
      // The inout MultiFab must have at least Periodic_BC_grow cells for the
      // algorithms taking care of periodic boundary conditions
    assert( inout.nGrow() >= LinOp_grow);
    
      // No coarsened boundary values, cannot apply inhomog at lev>0
    assert( !(level>0 && bc_mode == Inhomogeneous_BC) );
    
    int flagden = 1;        // fill in the bndry data and undrrelxr
    int flagbc  = 1;        // with values
    if (bc_mode == LinOp::Homogeneous_BC) flagbc = 0; // nodata if homog
    int nc = inout.nComp();

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
        const Array<BoundCond> &b = bgb.bndryConds(oitr());
        const Array<Real> &r = bgb.bndryLocs(oitr());
        FabSet &f = (*undrrelxr[level])[oitr()];
        int cdr(oitr());
        const FabSet& fs = bgb.bndryValues(oitr());
        //for(int gn = 0; gn < inout.length(); ++gn) {
        for(MultiFabIterator inoutmfi(inout); inoutmfi.isValid(); ++inoutmfi) {
          DependentFabSetIterator ffsi(inoutmfi, f);
          DependentFabSetIterator fsfsi(inoutmfi, fs);
            int gn = inoutmfi.index();
            assert(gbox[level]->get(inoutmfi.index()) == inoutmfi.validbox());
            const Mask& m = *maskvals[level][gn][oitr()];
            Real bcl(r[gn]);
            int bct(b[gn]);
            FORT_APPLYBC(
                &flagden, &flagbc, &maxorder,
                inoutmfi().dataPtr(), 
                ARLIM(inoutmfi().loVect()), ARLIM(inoutmfi().hiVect()),
                &cdr, &bct, &bcl,
                fsfsi().dataPtr(), 
                ARLIM(fsfsi().loVect()), ARLIM(fsfsi().hiVect()),
                m.dataPtr(), ARLIM(m.loVect()), ARLIM(m.hiVect()),
                ffsi().dataPtr(), ARLIM(ffsi().loVect()), ARLIM(ffsi().hiVect()),
                inoutmfi.validbox().loVect(), inoutmfi.validbox().hiVect(),
                &nc, h[level]);
        }
        oitr++;
    }
}
    
void
LinOp::residual(MultiFab &residL, const MultiFab &rhsL, MultiFab &solnL,
                int level, LinOp::BC_Mode bc_mode)
{
    apply(residL, solnL, level, bc_mode);
    //for(int gn = 0; gn < solnL.length(); ++gn) {
    for(MultiFabIterator solnLmfi(solnL); solnLmfi.isValid(); ++solnLmfi) {
      DependentMultiFabIterator residLmfi(solnLmfi, residL);
      DependentMultiFabIterator rhsLmfi(solnLmfi, rhsL);
        int nc = residL.nComp();
        assert(gbox[level]->get(solnLmfi.index()) == solnLmfi.validbox());
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
LinOp::smooth(MultiFab &solnL, const MultiFab &rhsL,
              int level, LinOp::BC_Mode bc_mode)
{
    for (int redBlackFlag = 0; redBlackFlag < 2; redBlackFlag++) {
        applyBC(solnL, level, bc_mode);
        Fsmooth(solnL, rhsL, level, redBlackFlag);
   }
}

Real
LinOp::norm(const MultiFab &in, int level) const
{
    Real norm = 0.0;
    //for(int gn = 0; gn < in.length(); ++gn) {
    for(ConstMultiFabIterator inmfi(in); inmfi.isValid(); ++inmfi) {
        int gn = inmfi.index();
        Real tnorm = inmfi().norm(gbox[level]->get(gn));
        norm += tnorm*tnorm;
    }
    ParallelDescriptor::ReduceRealSum(norm);
    return norm;
}

void
LinOp::clearToLevel(int level)
{
    int i;
    for(i = level+1; i < numLevels(); ++i) {
        delete undrrelxr[i].release();
        gbox[i] = 0;
    }
    h.resize(level+1);
    gbox.resize(level+1);
    undrrelxr.resize(level+1);
}

void
LinOp::computePeriodicIntersections(const Geometry&geo, const BoxArray& boxarr,
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
LinOp::prepareForLevel(int level)
{
    if(level == 0)
        return;
    LinOp::prepareForLevel(level-1);
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
    undrrelxr[level] = new BndryRegister(*gbox[level], 1, 0, 0, 1);
    maskvals.resize(level+1);
    maskvals[level].resize(gbox[level]->length());
    for(i=0; i < gbox[level]->length(); ++i) {
        maskvals[level][i].resize(2*BL_SPACEDIM, (Mask*)0);
            int m = 0;
            OrientationIter oitr;
            while ( oitr ) {
            Box bx_k = adjCell((*gbox[level])[i], oitr(), 1);
                maskvals[level][i][m] = new Mask(bx_k, 1);
                maskvals[level][i][m]->setVal(BndryData::not_covered);
            for(int gn = 0; gn < gbox[level]->length(); ++gn) {
                Box btmp = (*gbox[level])[gn];
                btmp &= bx_k;
                maskvals[level][i][m]->setVal(BndryData::covered, btmp,0);
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
                  Box btmp = (*gbox[level])[gn];
                  btmp &= curmask.box();
                  curmask.setVal(BndryData::covered, btmp,0);
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
LinOp::makeCoefficients(MultiFab& cs, const MultiFab &fn, int level)
{
    int nc = 1;
    const BoxArray& bxa = fn.boxArray();

      // Determine index type of incoming MultiFab
    IndexType itype = bxa[0].ixType();
    int cdir = -1;                // default to CELLTYPE
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

    int nComp=1;
    int nGrow=0;
    cs.define(d, nComp, nGrow, Fab_allocate);
    cs.setVal(0.0);

      // some abbreviations...
    BoxArray &mgbl = *(gbox[level]);
    //for(int gn = 0; gn < fn.length(); ++gn) {
    for(MultiFabIterator csmfi(cs); csmfi.isValid(); ++csmfi) {
      DependentMultiFabIterator fnmfi(csmfi, fn);
      //assert(mgbl[csmfi.index()] == csmfi.validbox());
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
