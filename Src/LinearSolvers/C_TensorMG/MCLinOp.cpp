
//
// $Id: MCLinOp.cpp,v 1.19 2004-01-07 21:18:43 car Exp $
//
// Differences from LinOp: den has nc components, bct has nc components.
//
#include <winstd.H>

#include <cstdlib>
#include <cmath>

#include <ParmParse.H>
#include <ParallelDescriptor.H>

#include "LO_BCTYPES.H"
#include "DivVis_F.H"
#include "MCLinOp.H"

bool MCLinOp::initialized = false;
int MCLinOp::def_harmavg  = 0;
int MCLinOp::def_verbose  = 0;
int MCLinOp::def_maxorder = 2;

#ifndef NDEBUG
//
// MCLinOp::applyBC fills MCLinOp_grow ghost cells with data expected in
// MCLinOp::apply() therefore, the incoming MultiFab to MCLinOp::applyBC()
// better have this many ghost allocated.
//
const int MCLinOp_grow = 1;
#endif

void
MCLinOp::initialize ()
{
    ParmParse pp("MCLp");
    pp.query("harmavg", def_harmavg);
    pp.query("v", def_verbose);
    if (ParallelDescriptor::IOProcessor() && def_verbose)
	std::cout << "def_harmavg = " << def_harmavg << '\n';
    initialized = true;
}

MCLinOp::MCLinOp (const BndryData& _bgb,
		  const Real       _h)
    :
    bgb(_bgb)
{
    BL_ASSERT (MCLinOp::bcComponentsNeeded() == bgb.nComp());
    Real __h[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        __h[i] = _h;
    }
    initConstruct(__h);
}

MCLinOp::MCLinOp (const BndryData& _bgb,
		  const Real*      _h)
    :
    bgb(_bgb)
{
    BL_ASSERT (MCLinOp::bcComponentsNeeded() == bgb.nComp());
    initConstruct(_h);
}

MCLinOp::~MCLinOp ()
{
    for (int i = 0; i < maskvals.size(); ++i)
    {
	for (int j=0; j < maskvals[i].size(); ++j)
        {
	    for (int k = 0; k < maskvals[i][j].size(); ++k) {
		delete maskvals[i][j][k];
	    }
	}
    }
}

MCLinOp::MCLinOp (const MCLinOp& _lp,
		  int            level)
    :
    bgb(_lp.bgb)
{
    BL_ASSERT(_lp.numLevels() > level);
    harmavg = _lp.harmavg;
    verbose = _lp.verbose;
    gbox.resize(1);
    gbox[0] = _lp.boxArray(level);
    geomarray.resize(1);
    geomarray[0] = bgb.getGeom();
    h.resize(1);
    h[0] = _lp.h[level];	// level should be prepared here.
    undrrelxr.resize(1);
    undrrelxr[0] = _lp.undrrelxr[level];
    tangderiv.resize(1);
    tangderiv[0] = _lp.tangderiv[level];
}

void
MCLinOp::initConstruct (const Real* _h)
{   
    if (!initialized)
	initialize();
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
    maskvals[level].resize(gbox[level].size());
    //
    // For each orientation, build NULL masks, then use distributed allocation
    //
    for (int i = 0; i < gbox[level].size(); i++)
    {
        maskvals[level][i].resize(2*BL_SPACEDIM, (Mask*)0);
    }

    const int myproc = ParallelDescriptor::MyProc();

    for (OrientationIter oitr; oitr; ++oitr)
    {
        Orientation face    = oitr();
        const FabSet& bndry = bgb[face];
        for (int i = 0; i < gbox[level].size(); i++)
        {
            if (bndry.DistributionMap()[i] == myproc)
            {
                const PArray<Mask>& pam = bgb.bndryMasks(face);
                BL_ASSERT(maskvals[level][i][face] == 0);
                maskvals[level][i][face] = new Mask(pam[i].box(), 1);
                maskvals[level][i][face]->copy(pam[i]);
            }
        }
    }
}

int
MCLinOp::bcComponentsNeeded()
{
  float d = (float) BL_SPACEDIM;
  int nDer = int(std::pow(d,d-1)+d);
  return nDer;
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

    inout.FillBoundary();
    prepareForLevel(level);

    geomarray[level].FillPeriodicBoundary(inout,0,nc);
    //
    // Fill boundary cells.
    //
    OrientationIter oitr;

    while (oitr)
    {
	const Array<Array<BoundCond> > &b = bgb.bndryConds(oitr());
	const Array<Real> &r = bgb.bndryLocs(oitr());
	FabSet& f  = (*undrrelxr[level])[oitr()];
	FabSet& td = (*tangderiv[level])[oitr()];
	int cdr(oitr());
	const FabSet& fs = bgb.bndryValues(oitr());
	int cdir = oitr().coordDir();
        for (MFIter inoutmfi(inout); inoutmfi.isValid(); ++inoutmfi)
        {
	    const int gn = inoutmfi.index();
	    BL_ASSERT(gbox[level][inoutmfi.index()] == inoutmfi.validbox());
	    Real bcl(r[gn]);
	    const int *bct = (const int *)b[gn].dataPtr();
	    const FArrayBox& fsfab = fs[inoutmfi.index()];
	    const Real* bcvalptr = fsfab.dataPtr();
            //
	    // Way external derivs stored.
            //
	    const Real* exttdptr = fsfab.dataPtr(numcomp); 
	    const int* fslo      = fsfab.loVect();
	    const int* fshi      = fsfab.hiVect();
	    FArrayBox& inoutfab  = inout[inoutmfi];
	    FArrayBox& denfab    = f[inoutmfi.index()];
	    FArrayBox& tdfab     = td[inoutmfi.index()];
#if BL_SPACEDIM==2
	    int perpdir;
	    if (cdir == 0)
                perpdir = 1;
	    else if (cdir == 1)
                perpdir = 0;
	    else
                BoxLib::Abort("MCLinOp::applyBC(): bad logic");

	    const Mask& m    = *maskvals[level][gn][oitr()];
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
		m.dataPtr(),    ARLIM(m.loVect()),    ARLIM(m.hiVect()),
		mphi.dataPtr(), ARLIM(mphi.loVect()), ARLIM(mphi.hiVect()),
		mplo.dataPtr(), ARLIM(mplo.loVect()), ARLIM(mplo.hiVect()),
		denfab.dataPtr(), 
		ARLIM(denfab.loVect()), ARLIM(denfab.hiVect()),
		exttdptr, ARLIM(fslo), ARLIM(fshi),
		tdfab.dataPtr(),ARLIM(tdfab.loVect()),ARLIM(tdfab.hiVect()),
		inoutmfi.validbox().loVect(), inoutmfi.validbox().hiVect(),
		&nc, h[level]);
#elif BL_SPACEDIM==3
	    const Mask& mn = *maskvals[level][gn][Orientation(1,Orientation::high)];
	    const Mask& me = *maskvals[level][gn][Orientation(0,Orientation::high)];
	    const Mask& mw = *maskvals[level][gn][Orientation(0,Orientation::low)];
	    const Mask& ms = *maskvals[level][gn][Orientation(1,Orientation::low)];
	    const Mask& mt = *maskvals[level][gn][Orientation(2,Orientation::high)];
	    const Mask& mb = *maskvals[level][gn][Orientation(2,Orientation::low)];
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
		inoutmfi.validbox().loVect(), inoutmfi.validbox().hiVect(),
		&nc, h[level]);
#endif
	}
	++oitr;
    }
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
	int nc = residL.nComp();
	FORT_RESIDL(
	    residL[solnLmfi].dataPtr(), 
            ARLIM(residL[solnLmfi].loVect()), ARLIM(residL[solnLmfi].hiVect()),
	    rhsL[solnLmfi].dataPtr(), 
            ARLIM(rhsL[solnLmfi].loVect()), ARLIM(rhsL[solnLmfi].hiVect()),
	    residL[solnLmfi].dataPtr(), 
            ARLIM(residL[solnLmfi].loVect()), ARLIM(residL[solnLmfi].hiVect()),
	    solnLmfi.validbox().loVect(), solnLmfi.validbox().hiVect(), &nc);
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
	delete undrrelxr[i].release();
	delete tangderiv[i].release();
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
    maskvals[level].resize(gbox[level].size());
    for (i = 0; i < gbox[level].size(); i++)
    {
        maskvals[level][i].resize(2*BL_SPACEDIM, (Mask*)0);
    }
    int myproc = ParallelDescriptor::MyProc();

    Array<IntVect> pshifts(27);

    for (OrientationIter oitr; oitr; ++oitr)
    {
        Orientation face = oitr();
        //
        // Use bgb's distribution map for masks.
        //
        for (FabSetIter bndryfsi(bgb[face]); bndryfsi.isValid(); ++bndryfsi)
	{
	    const int gn   = bndryfsi.index();
	    Box       bx_k = BoxLib::adjCell(gbox[level][gn], face, 1);
            //
	    // Extend box in directions orthogonal to face normal.
            //
	    for (int dir = 0; dir < BL_SPACEDIM; dir++)
            {
		if (dir == oitr())
                    continue;
		bx_k.grow(dir,1);
	    }
	    BL_ASSERT(maskvals[level][gn][face] == 0);
	    maskvals[level][gn][face] = new Mask(bx_k, 1);
	    Mask& curmask = *(maskvals[level][gn][face]);
	    curmask.setVal(BndryData::not_covered);
	    for (int gno = 0; gno < gbox[level].size(); ++gno)
            {
		Box btmp = gbox[level][gno] & bx_k;
		if (gno != gn  &&  btmp.ok())
		    curmask.setVal(BndryData::covered, btmp,0);
	    }
            //
	    // Now take care of periodic wraparounds.
            //
	    Geometry& curgeom = geomarray[level];
	    if (curgeom.isAnyPeriodic() && !curdomain.contains(bx_k))
	    {
		curgeom.periodicShift(curdomain, bx_k, pshifts);

		for (int iiv = 0; iiv<pshifts.size(); iiv++)
		{
		    curmask.shift(pshifts[iiv]);
		    for (int gno = 0; gno < gbox[level].size(); ++gno)
		    {
			Box btmp = gbox[level][gno] & curmask.box();
			curmask.setVal(BndryData::covered, btmp,0);
		    }

		    curmask.shift(-pshifts[iiv]);
		}
	    }
        }
    }
}

void
MCLinOp::makeCoefficients (MultiFab&       cs,
                           const MultiFab& fn,
                           int             level)
{
    int nc = fn.nComp();
    //
    // Determine index type of incoming MultiFab.
    //
    const IndexType iType(fn.boxArray()[0].ixType());
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
	switch(cdir)
        {
	case -1:
	    FORT_AVERAGECC(
		cs[csmfi].dataPtr(),
                ARLIM(cs[csmfi].loVect()), ARLIM(cs[csmfi].hiVect()),
		fn[csmfi].dataPtr(),
                ARLIM(fn[csmfi].loVect()), ARLIM(fn[csmfi].hiVect()),
		grids[csmfi.index()].loVect(),
                grids[csmfi.index()].hiVect(), &nc);
	    break;
	case 0:
	case 1:
	case 2:
	    if ( harmavg )
            {
		FORT_HARMONIC_AVERAGEEC(
		    cs[csmfi].dataPtr(), 
                    ARLIM(cs[csmfi].loVect()), ARLIM(cs[csmfi].hiVect()),
		    fn[csmfi].dataPtr(), 
                    ARLIM(fn[csmfi].loVect()), ARLIM(fn[csmfi].hiVect()),
		    grids[csmfi.index()].loVect(),
                    grids[csmfi.index()].hiVect(), &nc, &cdir);
	    }
            else
            {
		FORT_AVERAGEEC(
		    cs[csmfi].dataPtr(), 
                    ARLIM(cs[csmfi].loVect()), ARLIM(cs[csmfi].hiVect()),
		    fn[csmfi].dataPtr(), 
                    ARLIM(fn[csmfi].loVect()), ARLIM(fn[csmfi].hiVect()),
		    grids[csmfi.index()].loVect(),
                    grids[csmfi.index()].hiVect(), &nc, &cdir);
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
	os << "MCLinOp" << std::endl;
	os << "Grids: " << std::endl;
	for (int level = 0; level < lp.h.size(); ++level)
	{
	    os << " level = " << level << ": " << lp.gbox[level] << std::endl;
	}
	os << "Grid Spacing: " << std::endl;
	for (int level = 0; level < lp.h.size(); ++level)
	{
	    os << " level = " << level << ", dx = ";
	    for (int d =0; d < BL_SPACEDIM; ++d)
	    {
		os << lp.h[level][d] << "  ";
	    }
	    os << std::endl;
	}
	os << "Harmonic average? " << (lp.harmavg == 1 ? "yes" : "no") << std::endl;
	os << "Verbosity: " << lp.verbose << std::endl;
	os << "Max Order: " << lp.maxorder << std::endl;
    }

    if (ParallelDescriptor::IOProcessor())
    {
	os << "Masks:" << std::endl;
    }
    for (int level = 0; level < lp.h.size(); ++level)
    {
	if (ParallelDescriptor::IOProcessor())
	    os << "level = " << level << std::endl;

	for (int nproc = 0; nproc < ParallelDescriptor::NProcs(); ++nproc)
	{
	    if (nproc == ParallelDescriptor::MyProc())
	    {
		os << "Processor " << nproc << std::endl;
		for (OrientationIter oitr; oitr; ++oitr)
		{
		    Orientation face = oitr();
		    for (int i=0; i<lp.boxArray().size(); ++i)
		    {
			if (lp.maskvals[level][i][face])
			{
			    os << *lp.maskvals[level][i][face];
			}
		    }
		}
	    }
	}
    }    
    
    return os;
}

