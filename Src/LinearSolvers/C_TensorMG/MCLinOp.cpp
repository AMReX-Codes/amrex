#include <winstd.H>
#include <iostream>
#include <cstdlib>

#include <ParmParse.H>
#include <ParallelDescriptor.H>

#include <LO_BCTYPES.H>
#include <MCLO_F.H>
#include <MCLinOp.H>

namespace
{
    bool initialized = false;
}
//
// Set default values for these in Initialize()!!!
//
int MCLinOp::def_harmavg;
int MCLinOp::def_verbose;
int MCLinOp::def_maxorder;
int MCLinOp::def_ncomp = BL_SPACEDIM;

#ifndef NDEBUG
//
// MCLinOp::applyBC fills MCLinOp_grow ghost cells with data expected in
// MCLinOp::apply() therefore, the incoming MultiFab to MCLinOp::applyBC()
// better have this many ghost allocated.
//
const int MCLinOp_grow = 1;
#endif

void
MCLinOp::Initialize ()
{
    if (initialized) return;
    MCLinOp::def_harmavg  = 0;
    MCLinOp::def_verbose  = 0;
    MCLinOp::def_maxorder = 2;

    ParmParse pp("MCLp");

    pp.query("harmavg", def_harmavg);
    pp.query("v",       def_verbose);
    pp.query("maxorder",def_maxorder);

    if (ParallelDescriptor::IOProcessor() && def_verbose)
	std::cout << "def_harmavg = " << def_harmavg << '\n';

    BoxLib::ExecOnFinalize(MCLinOp::Finalize);

    initialized = true;
}

void
MCLinOp::Finalize ()
{
    initialized = false;
}

MCLinOp::MCLinOp (const BndryData& _bgb,
		  const Real       _h,
		  int              _nc)
    : numcomp(_nc), bgb(_bgb)
{
    BL_ASSERT (MCLinOp::bcComponentsNeeded(numcomp) == bgb.nComp());
    Real __h[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        __h[i] = _h;
    }
    initConstruct(__h);
}

MCLinOp::MCLinOp (const BndryData& _bgb,
		  const Real*      _h,
		  int              _nc)
    : numcomp(_nc), bgb(_bgb)
{
    BL_ASSERT (MCLinOp::bcComponentsNeeded(numcomp) == bgb.nComp());
    initConstruct(_h);
}

MCLinOp::~MCLinOp ()
{
    for (int i = 0, N = maskvals.size(); i < N; ++i)
    {
        for (std::map<int,MaskTuple>::iterator it = maskvals[i].begin(),
                 End = maskvals[i].end();
             it != End;
             ++it)
        {
            MaskTuple& a = it->second;
            for (int k = 0; k < 2*BL_SPACEDIM; ++k)
                delete a[k];
        }
    }
}

void
MCLinOp::initConstruct (const Real* _h)
{   
    Initialize();
    //
    // We'll reserve() space to cut down on copying during resize()s.
    //
    const int N = 10;

    h.reserve(N);
    gbox.reserve(N);
    undrrelxr.reserve(N);
    tangderiv.reserve(N);
    maskvals.reserve(N);
    geomarray.reserve(N);

    harmavg = def_harmavg;
    verbose = def_verbose;
    gbox.resize(1);
    int level = 0;
    gbox[level] = bgb.boxes();
    geomarray.resize(1);
    geomarray[level] = bgb.getGeom();
    h.resize(1);
    maxorder = def_maxorder;
    for (int i = 0; i < BL_SPACEDIM; ++i)
    {
	h[level][i] = _h[i];
    }
    maskvals.resize(1);
    //
    // For each orientation, build NULL masks, then use distributed allocation.
    // We note that all orientations of the FabSets have the same distribution.
    // We'll use the low 0 side as the model.
    //
    for (FabSetIter bndryfsi(bgb[Orientation(0,Orientation::low)]);
         bndryfsi.isValid();
         ++bndryfsi)
    {
        const int        i   = bndryfsi.index();
        MaskTuple&       ma  =  maskvals[level][i];
        const MaskTuple& bdm = bgb.bndryMasks(i);

        for (OrientationIter oitr; oitr; ++oitr)
        {
            const Orientation face  = oitr();
            const Mask*       m     = bdm[face];
            ma[face] = new Mask(m->box(),1);
            ma[face]->copy(*m);
        }
    }
}

int
MCLinOp::bcComponentsNeeded(int ncomp)
{
  int nc;
#if (BL_SPACEDIM==2)
  nc = ncomp * 2; // Tangential derivatives for each comp
#else
  nc = ncomp * (1 + BL_SPACEDIM); // D-1 tang derivs for ea, but waste one slot to simplify indexing
#endif
  return nc;
}

void
MCLinOp::apply (MultiFab& out,
		MultiFab& in,
		int       level,
		MCBC_Mode bc_mode)
{
    applyBC(in,level,bc_mode);
    Fapply(out,in,level);
}

void
MCLinOp::applyBC (MultiFab& inout,
		  int       level,
		  MCBC_Mode bc_mode)
{
    //
    // The inout MultiFab must have at least MCLinOp_grow ghost cells
    // for applyBC()
    //
    BL_ASSERT(inout.nGrow() >= MCLinOp_grow);
    //
    // The inout MultiFab must have at least Periodic_BC_grow cells for the
    // algorithms taking care of periodic boundary conditions.
    //
    BL_ASSERT(inout.nGrow() >= MCLinOp_grow);
    //
    // No coarsened boundary values, cannot apply inhomog at lev>0.
    //
    BL_ASSERT(!(level>0 && bc_mode == MCInhomogeneous_BC));
    
    int flagden = 1;	// fill in the bndry data and undrrelxr
    int flagbc  = 1;	// with values
    if (bc_mode == MCHomogeneous_BC)
        flagbc = 0; // nodata if homog
    int nc = inout.nComp();
    BL_ASSERT(nc == numcomp );

    inout.setBndry(-1.e30);
    inout.FillBoundary();
    prepareForLevel(level);

    geomarray[level].FillPeriodicBoundary(inout,0,nc);
    //
    // Fill boundary cells.
    //
    const int N = inout.IndexMap().size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
    {
        const int gn = inout.IndexMap()[i];

        BL_ASSERT(gbox[level][gn] == inout.box(gn));

        const BndryData::RealTuple&      bdl = bgb.bndryLocs(gn);
        const Array< Array<BoundCond> >& bdc = bgb.bndryConds(gn);
        const MaskTuple&                 msk = maskvals[level][gn];

        for (OrientationIter oitr; oitr; ++oitr)
        {
            const Orientation face = oitr();
            FabSet& f  = (*undrrelxr[level])[face];
            FabSet& td = (*tangderiv[level])[face];
            int cdr(face);
            const FabSet& fs = bgb.bndryValues(face);
	    Real bcl = bdl[face];
            const Array<BoundCond>& bc = bdc[face];
	    const int *bct = (const int*) bc.dataPtr();
	    const FArrayBox& fsfab = fs[gn];
	    const Real* bcvalptr = fsfab.dataPtr();
            //
	    // Way external derivs stored.
            //
	    const Real* exttdptr = fsfab.dataPtr(numcomp); 
	    const int* fslo      = fsfab.loVect();
	    const int* fshi      = fsfab.hiVect();
	    FArrayBox& inoutfab  = inout[gn];
	    FArrayBox& denfab    = f[gn];
	    FArrayBox& tdfab     = td[gn];
#if BL_SPACEDIM==2
            int cdir = face.coordDir(), perpdir = -1;
	    if (cdir == 0)
                perpdir = 1;
	    else if (cdir == 1)
                perpdir = 0;
	    else
                BoxLib::Abort("MCLinOp::applyBC(): bad logic");

	    const Mask& m    = *msk[face];
	    const Mask& mphi = *msk[Orientation(perpdir,Orientation::high)];
	    const Mask& mplo = *msk[Orientation(perpdir,Orientation::low)];
	    FORT_APPLYBC(
		&flagden, &flagbc, &maxorder,
		inoutfab.dataPtr(), 
                ARLIM(inoutfab.loVect()), ARLIM(inoutfab.hiVect()),
		&cdr, bct, &bcl,
		bcvalptr, ARLIM(fslo), ARLIM(fshi),
		m.dataPtr(),    ARLIM(m.loVect()),    ARLIM(m.hiVect()),
		mphi.dataPtr(), ARLIM(mphi.loVect()), ARLIM(mphi.hiVect()),
		mplo.dataPtr(), ARLIM(mplo.loVect()), ARLIM(mplo.hiVect()),
		denfab.dataPtr(), 
		ARLIM(denfab.loVect()), ARLIM(denfab.hiVect()),
		exttdptr, ARLIM(fslo), ARLIM(fshi),
		tdfab.dataPtr(),ARLIM(tdfab.loVect()),ARLIM(tdfab.hiVect()),
		inout.box(gn).loVect(), inout.box(gn).hiVect(),
		&nc, h[level]);
#elif BL_SPACEDIM==3
	    const Mask& mn = *msk[Orientation(1,Orientation::high)];
	    const Mask& me = *msk[Orientation(0,Orientation::high)];
	    const Mask& mw = *msk[Orientation(0,Orientation::low)];
	    const Mask& ms = *msk[Orientation(1,Orientation::low)];
	    const Mask& mt = *msk[Orientation(2,Orientation::high)];
	    const Mask& mb = *msk[Orientation(2,Orientation::low)];
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
		inout.box(gn).loVect(), inout.box(gn).hiVect(),
		&nc, h[level]);
#endif
	}
    }

#if 0
  // This "probably" works, but is not strictly needed just because of the way Bill
  // coded up the tangential derivative stuff.  It's handy code though, so I want to
  // keep it around/

  // Clean up corners:
  // The problem here is that APPLYBC fills only grow cells normal to the boundary.
  // As a result, any corner cell on the boundary (either coarse-fine or fine-fine)
  // is not filled.  For coarse-fine, the operator adjusts itself, sliding away from
  // the box edge to avoid referencing that corner point.  On the physical boundary
  // though, the corner point is needed.  Particularly if a fine-fine boundary intersects
  // the physical boundary, since we want the stencil to be independent of the box
  // blocking.  FillBoundary operations wont fix the problem because the "good"
  // data we need is living in the grow region of adjacent fabs.  So, here we play
  // the usual games to treat the newly filled grow cells as "valid" data.

  // Note that we only need to do something where the grids touch the physical boundary.

  const Geometry& geomlev = geomarray[level];
  const BoxArray& grids = inout.boxArray();
  const Box& domain = geomlev.Domain();
  int nGrow = 1;
  int src_comp = 0;
  int num_comp = BL_SPACEDIM;


  // Lets do a quick check to see if we need to do anything at all here
  BoxArray BIGba = BoxArray(grids).grow(nGrow);

  if (! (domain.contains(BIGba.minimalBox())) ) {

    BoxArray boundary_pieces;
    Array<int> proc_idxs;
    Array<Array<int> > old_to_new(grids.size());
    const DistributionMapping& dmap=inout.DistributionMap();

    for (int d=0; d<BL_SPACEDIM; ++d) {
      if (! (geomlev.isPeriodic(d)) ) {

        BoxArray gba = BoxArray(grids).grow(d,nGrow);
        for (int i=0; i<gba.size(); ++i) {
          BoxArray new_pieces = BoxLib::boxComplement(gba[i],domain);
          int size_new = new_pieces.size();
          if (size_new>0) {
            int size_old = boundary_pieces.size();
            boundary_pieces.resize(size_old+size_new);
            proc_idxs.resize(boundary_pieces.size());
            for (int j=0; j<size_new; ++j) {
              boundary_pieces.set(size_old+j,new_pieces[j]);
              proc_idxs[size_old+j] = dmap[i];
              old_to_new[i].push_back(size_old+j);
            }
          }
        }
      }
    }

    proc_idxs.push_back(ParallelDescriptor::MyProc());

    MultiFab boundary_data(boundary_pieces,num_comp,nGrow,
                           DistributionMapping(proc_idxs));

    for (MFIter mfi(inout); mfi.isValid(); ++mfi) {
      const FArrayBox& src_fab = inout[mfi];
      for (int j=0; j<old_to_new[mfi.index()].size(); ++j) {
        int new_box_idx = old_to_new[mfi.index()][j];
        boundary_data[new_box_idx].copy(src_fab,src_comp,0,num_comp);
      }
    }

    boundary_data.FillBoundary();

    // Use a hacked Geometry object to handle the periodic intersections for us.
    // Here, the "domain" is the plane of cells on non-periodic boundary faces.
    // and there may be cells over the periodic boundary in the remaining directions.
    // We do a Geometry::PFB on each non-periodic face to sync these up.
    if (geomlev.isAnyPeriodic()) {
      Array<int> is_per(BL_SPACEDIM,0);
      for (int d=0; d<BL_SPACEDIM; ++d) {
        is_per[d] = geomlev.isPeriodic(d);
      }
      for (int d=0; d<BL_SPACEDIM; ++d) {
        if (! is_per[d]) {
          Box tmpLo = BoxLib::adjCellLo(geomlev.Domain(),d,1);
          Geometry tmpGeomLo(tmpLo,&(geomlev.ProbDomain()),(int)geomlev.Coord(),is_per.dataPtr());
          tmpGeomLo.FillPeriodicBoundary(boundary_data);

          Box tmpHi = BoxLib::adjCellHi(geomlev.Domain(),d,1);
          Geometry tmpGeomHi(tmpHi,&(geomlev.ProbDomain()),(int)geomlev.Coord(),is_per.dataPtr());
          tmpGeomHi.FillPeriodicBoundary(boundary_data);
        }
      }
    }

    for (MFIter mfi(inout); mfi.isValid(); ++mfi) {
      int idx = mfi.index();
      FArrayBox& dst_fab = inout[mfi];
      for (int j=0; j<old_to_new[idx].size(); ++j) {
        int new_box_idx = old_to_new[mfi.index()][j];
        const FArrayBox& src_fab = boundary_data[new_box_idx];
        const Box& src_box = src_fab.box();

        BoxArray pieces_outside_domain = BoxLib::boxComplement(src_box,domain);
        for (int k=0; k<pieces_outside_domain.size(); ++k) {
          const Box& outside = pieces_outside_domain[k] & dst_fab.box();
          if (outside.ok()) {
            dst_fab.copy(src_fab,outside,0,outside,src_comp,num_comp);
          }
        }
      }
    }
  }
#endif
}
    
void
MCLinOp::residual (MultiFab&       residL,
		   const MultiFab& rhsL,
		   MultiFab&       solnL,
		   int             level,
		   MCBC_Mode       bc_mode)
{
    apply(residL, solnL, level, bc_mode);

    for (MFIter solnLmfi(solnL); solnLmfi.isValid(); ++solnLmfi)
    {
	int              nc     = residL.nComp();
        const Box&       vbox   = solnLmfi.validbox();
        FArrayBox&       resfab = residL[solnLmfi];
        const FArrayBox& rhsfab = rhsL[solnLmfi];
	FORT_RESIDL(
	    resfab.dataPtr(), 
            ARLIM(resfab.loVect()), ARLIM(resfab.hiVect()),
	    rhsfab.dataPtr(), 
            ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
	    resfab.dataPtr(), 
            ARLIM(resfab.loVect()), ARLIM(resfab.hiVect()),
	    vbox.loVect(), vbox.hiVect(), &nc);
    }
}

void
MCLinOp::smooth (MultiFab&       solnL,
		 const MultiFab& rhsL,
		 int             level,
		 MCBC_Mode       bc_mode)
{
    for (int phaseflag = 0; phaseflag < numphase; phaseflag++)
    {
	applyBC(solnL, level, bc_mode);
	Fsmooth(solnL, rhsL, level, phaseflag);
    }
}

Real
MCLinOp::norm (const MultiFab& in,
	       int             level) const
{
    Real norm = 0.0;
    for (MFIter inmfi(in); inmfi.isValid(); ++inmfi)
    {
        Real tnorm = in[inmfi].norm(gbox[level][inmfi.index()]);
	norm += tnorm*tnorm;
    }
    ParallelDescriptor::ReduceRealSum(norm);
    return norm;
}

void
MCLinOp::clearToLevel (int level)
{
    for (int i = level+1; i < numLevels(); ++i)
    {
	gbox[i].clear();
    }
    h.resize(level+1);
    gbox.resize(level+1);
    undrrelxr.resize(level+1);
    tangderiv.resize(level+1);
}

void
MCLinOp::prepareForLevel (int level)
{
    if (level == 0) return;

    MCLinOp::prepareForLevel(level-1);

    if (h.size() > level) return;
    //
    // Assume from here down that this is a new level one coarser than existing
    //
    BL_ASSERT(h.size() == level);
    h.resize(level+1);
    int i;
    for (i = 0; i < BL_SPACEDIM; ++i)
	h[level][i] = h[level-1][i]*2.0;

    geomarray.resize(level+1);
    Box curdomain = Box( geomarray[level-1].Domain() ).coarsen(2);
    geomarray[level].define( curdomain );
    //
    // Add a box to the new coarser level (assign removes old BoxArray)
    //
    gbox.resize(level+1);
    gbox[level] = BoxArray(gbox[level-1]).coarsen(2);
    //
    // Add the BndryRegister of relax values to the new coarser level.
    //
    BL_ASSERT(undrrelxr.size() == level);
    undrrelxr.resize(level+1);
    undrrelxr[level] = new BndryRegister(gbox[level], 1, 0, 0, numcomp);
    //
    // Add the BndryRegister to hold tagential derivatives to the new
    // coarser level.
    //
    BL_ASSERT(tangderiv.size() == level);
    tangderiv.resize(level+1);
    //
    // Figure out how many components.
    //
    const FabSet& samplefs = (*tangderiv[level-1])[Orientation(0,Orientation::low)];
    tangderiv[level] = new BndryRegister(gbox[level],0,1,0,samplefs.nComp());
    //
    // Add an Array of Array of maskvals to the new coarser level
    // For each orientation, build NULL masks, then use distributed allocation
    // Initial masks for coarse levels, ignore outside_domain possibility since
    // we always solve homogeneous equation on coarse levels.
    //
    BL_ASSERT(maskvals.size() == level);
    maskvals.resize(level+1);

    Array<IntVect> pshifts(27);

    std::vector< std::pair<int,Box> > isects;
    //
    // Use bgb's distribution map for masks.
    // We note that all orientations of the FabSets have the same distribution.
    // We'll use the low 0 side as the model.
    //
    for (FabSetIter bndryfsi(bgb[Orientation(0,Orientation::low)]);
         bndryfsi.isValid();
         ++bndryfsi)
    {
        const int  gn  = bndryfsi.index();
        MaskTuple& msk = maskvals[level][gn];

        for (OrientationIter oitr; oitr; ++oitr)
        {
            const Orientation face = oitr();
	    Box               bx_k = BoxLib::adjCell(gbox[level][gn], face, 1);
            //
	    // Extend box in directions orthogonal to face normal.
            //
	    for (int dir = 0; dir < BL_SPACEDIM; dir++)
		if (dir != face)
                    bx_k.grow(dir,1);

	    msk[face] = new Mask(bx_k, 1);
	    Mask& curmask = *(msk[face]);
	    curmask.setVal(BndryData::not_covered);

            gbox[level].intersections(bx_k,isects);
            for (int ii = 0, N = isects.size(); ii < N; ii++)
                if (isects[ii].first != gn)
		    curmask.setVal(BndryData::covered, isects[ii].second, 0);
            //
	    // Now take care of periodic wraparounds.
            //
	    Geometry& curgeom = geomarray[level];

	    if (curgeom.isAnyPeriodic() && !curdomain.contains(bx_k))
	    {
		curgeom.periodicShift(curdomain, bx_k, pshifts);

                for (int iiv = 0, M = pshifts.size(); iiv < M; iiv++)
		{
		    curmask.shift(pshifts[iiv]);
                    gbox[level].intersections(curmask.box(),isects);
                    for (int ii = 0, N = isects.size(); ii < N; ii++)
                         curmask.setVal(BndryData::covered, isects[ii].second, 0);
		    curmask.shift(-pshifts[iiv]);
		}
	    }
        }
    }

    gbox[level].clear_hash_bin();
}

void
MCLinOp::makeCoefficients (MultiFab&       cs,
                           const MultiFab& fn,
                           int             level)
{
    const int nc = fn.nComp();
    //
    // Determine index type of incoming MultiFab.
    //
    const IndexType iType(fn.boxArray().ixType());
    const IndexType cType(D_DECL(IndexType::CELL, IndexType::CELL, IndexType::CELL));
    const IndexType xType(D_DECL(IndexType::NODE, IndexType::CELL, IndexType::CELL));
    const IndexType yType(D_DECL(IndexType::CELL, IndexType::NODE, IndexType::CELL));
#if (BL_SPACEDIM == 3)    
    const IndexType zType(D_DECL(IndexType::CELL, IndexType::CELL, IndexType::NODE));
#endif
    int cdir;
    if (iType == cType)
    {
        cdir = -1;
    }
    else if (iType == xType)
    {
        cdir = 0;
    }
    else if (iType == yType)
    {
        cdir = 1;
    }
#if (BL_SPACEDIM == 3)
    else if (iType == zType)
    {
        cdir = 2;
    }
#endif
    else
        BoxLib::Abort("MCLinOp::makeCoeffients(): Bad index type");
    
    BoxArray d(gbox[level]);
    if (cdir >= 0)
	d.surroundingNodes(cdir);

    int nGrow=0;
    cs.define(d, nc, nGrow, Fab_allocate);
    cs.setVal(0.0);

    const BoxArray& grids = gbox[level];

    for (MFIter csmfi(cs); csmfi.isValid(); ++csmfi)
    {
        const Box&       grd   = grids[csmfi.index()];
        FArrayBox&       csfab = cs[csmfi];
        const FArrayBox& fnfab = fn[csmfi];

	switch(cdir)
        {
	case -1:
	    FORT_AVERAGECC(
		csfab.dataPtr(),
                ARLIM(csfab.loVect()), ARLIM(csfab.hiVect()),
		fnfab.dataPtr(),
                ARLIM(fnfab.loVect()), ARLIM(fnfab.hiVect()),
		grd.loVect(),
                grd.hiVect(), &nc);
	    break;
	case 0:
	case 1:
	case 2:
	    if ( harmavg )
            {
		FORT_HARMONIC_AVERAGEEC(
		    csfab.dataPtr(), 
                    ARLIM(csfab.loVect()), ARLIM(csfab.hiVect()),
		    fnfab.dataPtr(), 
                    ARLIM(fnfab.loVect()), ARLIM(fnfab.hiVect()),
		    grd.loVect(),
                    grd.hiVect(), &nc, &cdir);
	    }
            else
            {
		FORT_AVERAGEEC(
		    csfab.dataPtr(), 
                    ARLIM(csfab.loVect()), ARLIM(csfab.hiVect()),
		    fnfab.dataPtr(), 
                    ARLIM(fnfab.loVect()), ARLIM(fnfab.hiVect()),
		    grd.loVect(),
                    grd.hiVect(), &nc, &cdir);
	    }
	    break;
	default:
	    BoxLib::Error("MCLinOp::makeCoeffients(): bad coefficient coarsening direction!");
	}
    }
}

std::ostream&
operator<< (std::ostream&  os,
            const MCLinOp& lp)
{
    if (ParallelDescriptor::IOProcessor())
    {
	os << "MCLinOp" << '\n';
	os << "Grids: " << '\n';
	for (int level = 0; level < lp.h.size(); ++level)
	{
	    os << " level = " << level << ": " << lp.gbox[level] << '\n';
	}
	os << "Grid Spacing: " << '\n';
	for (int level = 0; level < lp.h.size(); ++level)
	{
	    os << " level = " << level << ", dx = ";
	    for (int d =0; d < BL_SPACEDIM; ++d)
	    {
		os << lp.h[level][d] << "  ";
	    }
	    os << '\n';
	}
	os << "Harmonic average? " << (lp.harmavg == 1 ? "yes" : "no") << '\n';
	os << "Verbosity: " << lp.verbose << '\n';
	os << "Max Order: " << lp.maxorder << '\n';
    }

    if (ParallelDescriptor::IOProcessor())
	os << "Masks:" << '\n';

    for (int level = 0; level < lp.h.size(); ++level)
    {
	if (ParallelDescriptor::IOProcessor())
	    os << "level = " << level << '\n';

        const std::map<int,MCLinOp::MaskTuple>& m = lp.maskvals[level];

	for (int nproc = 0; nproc < ParallelDescriptor::NProcs(); ++nproc)
	{
	    if (nproc == ParallelDescriptor::MyProc())
	    {
		os << "Processor " << nproc << '\n';

		for (OrientationIter oitr; oitr; ++oitr)
		{
		    const Orientation face = oitr();

                    for (std::map<int,MCLinOp::MaskTuple>::const_iterator it = m.begin(),
                             End = m.end();
                         it != End;
                         ++it)
                    {
                        os << *(it->second[face]);
                    }
		}
	    }
	}
    }    
    
    return os;
}

