#include <iostream>
#include <cstdlib>

#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>

#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MCLO_F.H>
#include <AMReX_MCLinOp.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

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

//
// MCLinOp::applyBC fills MCLinOp_grow ghost cells with data expected in
// MCLinOp::apply() therefore, the incoming MultiFab to MCLinOp::applyBC()
// better have this many ghost allocated.
//
namespace {
    const int MCLinOp_grow = 1;
}

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

    if (ParallelDescriptor::IOProcessor() && def_verbose) {
        amrex::Print() << "def_harmavg = " << def_harmavg << '\n';
    }

    amrex::ExecOnFinalize(MCLinOp::Finalize);

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
    Real _hh[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        _hh[i] = _h;
    }
    initConstruct(_hh);
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
    maskvals[0].resize(2*BL_SPACEDIM);

    for (OrientationIter oitr; oitr; ++oitr)
    {
	const Orientation face = oitr();
	const MultiMask& m = bgb.bndryMasks(face);
	maskvals[0][face].define(m.boxArray(), m.DistributionMap(), 1);
	MultiMask::Copy(maskvals[0][face], m);
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

    prepareForLevel(level);

    inout.FillBoundary(geomarray[level].periodicity());

    //
    // Fill boundary cells.
    //
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(inout); mfi.isValid(); ++mfi)
    {
        const int gn = mfi.index();

	const Box& iobx = inout.box(gn);

        BL_ASSERT(gbox[level][gn] == inout.box(gn));

        const BndryData::RealTuple&      bdl = bgb.bndryLocs(gn);
        const Vector< Vector<BoundCond> >& bdc = bgb.bndryConds(gn);

        for (OrientationIter oitr; oitr; ++oitr)
        {
            const Orientation face = oitr();
            FabSet& f  = undrrelxr[level][face];
            FabSet& td = tangderiv[level][face];
            int cdr(face);
            const FabSet& fs = bgb.bndryValues(face);
	    Real bcl = bdl[face];
            const Vector<BoundCond>& bc = bdc[face];
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
                amrex::Abort("MCLinOp::applyBC(): bad logic");

            const Mask& m    = maskvals[level][face][mfi];
	    const Mask& mphi = maskvals[level][Orientation(perpdir,Orientation::high)][mfi];
	    const Mask& mplo = maskvals[level][Orientation(perpdir,Orientation::low) ][mfi];
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
		iobx.loVect(), iobx.hiVect(),
		&nc, h[level].data());
#elif BL_SPACEDIM==3
	    const Mask& mn = maskvals[level][Orientation(1,Orientation::high)][mfi];
	    const Mask& me = maskvals[level][Orientation(0,Orientation::high)][mfi];
	    const Mask& mw = maskvals[level][Orientation(0,Orientation::low) ][mfi];
	    const Mask& ms = maskvals[level][Orientation(1,Orientation::low) ][mfi];
	    const Mask& mt = maskvals[level][Orientation(2,Orientation::high)][mfi];
	    const Mask& mb = maskvals[level][Orientation(2,Orientation::low) ][mfi];
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
		iobx.loVect(), iobx.hiVect(),
		&nc, h[level].data());
#endif
	}
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
    MultiFab::Xpay(residL, -1.0, rhsL, 0, 0, residL.nComp(), 0);
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
    Real nm = MultiFab::Dot(in, 0, in, 0, in.nComp(), 0);
    return nm;
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
    undrrelxr[level].define(gbox[level], DistributionMap(), 1, 0, 0, numcomp);
    //
    // Add the BndryRegister to hold tagential derivatives to the new
    // coarser level.
    //
    BL_ASSERT(tangderiv.size() == level);
    tangderiv.resize(level+1);
    //
    // Figure out how many components.
    //
    const FabSet& samplefs = tangderiv[level-1][Orientation(0,Orientation::low)];
    tangderiv[level].define(gbox[level], DistributionMap(), 0,1,0,samplefs.nComp());
    //
    // Add an Array of Array of maskvals to the new coarser level
    // For each orientation, build NULL masks, then use distributed allocation
    // Initial masks for coarse levels, ignore outside_domain possibility since
    // we always solve homogeneous equation on coarse levels.
    //
    BL_ASSERT(maskvals.size() == level);
    maskvals.resize(level+1);
    maskvals[level].resize(2*BL_SPACEDIM);

    for (OrientationIter fi; fi; ++fi)
    {
        Orientation face = fi();
	maskvals[level][face].define(gbox[level], 
				     DistributionMap(),
				     geomarray[level],
				     face, 0, 1, 1, 1, true);
    }
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
    const IndexType cType(AMREX_D_DECL(IndexType::CELL, IndexType::CELL, IndexType::CELL));
    const IndexType xType(AMREX_D_DECL(IndexType::NODE, IndexType::CELL, IndexType::CELL));
    const IndexType yType(AMREX_D_DECL(IndexType::CELL, IndexType::NODE, IndexType::CELL));
#if (BL_SPACEDIM == 3)    
    const IndexType zType(AMREX_D_DECL(IndexType::CELL, IndexType::CELL, IndexType::NODE));
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
        amrex::Abort("MCLinOp::makeCoeffients(): Bad index type");
    
    BoxArray d(gbox[level]);
    if (cdir >= 0)
	d.surroundingNodes(cdir);

    int nGrow=0;
    cs.define(d, DistributionMap(), nc, nGrow);
    cs.setVal(0.0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter csmfi(cs,true); csmfi.isValid(); ++csmfi)
    {
        const Box&       bx    = csmfi.tilebox();
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
		bx.loVect(),
                bx.hiVect(), &nc);
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
		    bx.loVect(),
                    bx.hiVect(), &nc, &cdir);
	    }
            else
            {
		FORT_AVERAGEEC(
		    csfab.dataPtr(), 
                    ARLIM(csfab.loVect()), ARLIM(csfab.hiVect()),
		    fnfab.dataPtr(), 
                    ARLIM(fnfab.loVect()), ARLIM(fnfab.hiVect()),
		    bx.loVect(),
                    bx.hiVect(), &nc, &cdir);
	    }
	    break;
	default:
	    amrex::Error("MCLinOp::makeCoeffients(): bad coefficient coarsening direction!");
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

	for (int nproc = 0; nproc < ParallelDescriptor::NProcs(); ++nproc)
	{
	    if (nproc == ParallelDescriptor::MyProc())
	    {
		os << "Processor " << nproc << '\n';

		for (OrientationIter oitr; oitr; ++oitr)
		{
		    const Orientation face = oitr();

		    for (MultiMaskIter mmi(lp.maskvals[level][face]); mmi.isValid(); ++mmi)
		    {
                        os << lp.maskvals[level][face][mmi];
                    }
		}
	    }
	}
    }    
    
    return os;
}

}
