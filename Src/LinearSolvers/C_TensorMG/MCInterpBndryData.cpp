#ifdef BL_USE_NEW_HFILES
#include <cmath>
#else
#include <math.h>
#endif

#include <LO_BCTYPES.H>
#include <MCInterpBndryData.H>
#include <MCINTERPBNDRYDATA_F.H>

static BDInterpFunc* bdfunc[2*BL_SPACEDIM];
static BDPhysDerivative* bdider[2*BL_SPACEDIM];
static int bdfunc_set = 0;

static void bdfunc_init()
{
    Orientation xloface(0,Orientation::low);
    Orientation xhiface(0,Orientation::high);
    bdfunc[xloface] = FORT_BDINTERPXLO;
    bdfunc[xhiface] = FORT_BDINTERPXHI;
    bdider[xloface] = FORT_BDIDERIVXLO;
    bdider[xhiface] = FORT_BDIDERIVXHI;
#if (BL_SPACEDIM > 1)
    Orientation yloface(1,Orientation::low);
    Orientation yhiface(1,Orientation::high);
    bdfunc[yloface] = FORT_BDINTERPYLO;
    bdfunc[yhiface] = FORT_BDINTERPYHI;
    bdider[yloface] = FORT_BDIDERIVYLO;
    bdider[yhiface] = FORT_BDIDERIVYHI;
#endif
#if (BL_SPACEDIM > 2)
    Orientation zloface(2,Orientation::low);
    Orientation zhiface(2,Orientation::high);
    bdfunc[zloface] = FORT_BDINTERPZLO;
    bdfunc[zhiface] = FORT_BDINTERPZHI;
    bdider[zloface] = FORT_BDIDERIVZLO;
    bdider[zhiface] = FORT_BDIDERIVZHI;
#endif
}

#if (BL_SPACEDIM == 2)
#define NUMDERIV 2
#endif

#if (BL_SPACEDIM == 3)
#define NUMDERIV 9
#endif

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();
#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

MCInterpBndryData::MCInterpBndryData (const BoxArray& _grids,
				      int             _ncomp,
				      const Geometry& geom)
    : BndryData(_grids,_ncomp,geom)
{}

// At the coarsest level the bndry values are taken from adjacent grids.
void
MCInterpBndryData::setBndryValues(const MultiFab&     mf,
				  int                 mf_start,
				  int                 bnd_start,
				  int                 num_comp,
				  const Array<BCRec>& bc )
{
    if (! bdfunc_set) bdfunc_init();
    
      // check that boxarrays are identical
    assert( grids.ready() );
    assert( grids == mf.boxArray() );

      // set bndry flags and locations
    int ref_ratio = 1;
    setBndryConds(bc, geom, ref_ratio);

    // find the dx's
    const Real* h = geom.CellSize();

    // HACK: cast away const to satisfy incomplete BoxLib interface
    for(MultiFabIterator mfi((MultiFab&)mf); mfi.isValid(); ++mfi)
    {
	assert(grids[mfi.index()] == mfi.validbox());
        const Box& bx = mfi.validbox();
        for(OrientationIter fi; fi; ++fi) {
            Orientation face(fi());
	    int dir = face.coordDir();

	    // physical bndry, copy from grid
            if(bx[face] == geom.Domain()[face])
	    {
		// load up hfine with perpindicular h's
		Real hfine[BL_SPACEDIM];
		int kdir = 0;
		for(int idir=0; idir<BL_SPACEDIM; ++idir){
		    if( idir == dir ) continue;
		    hfine[kdir++] = h[idir];
		}

		DependentFabSetIterator bfsi(mfi,bndry[face]);
		
		// copy and compute deriv
		bdider[face](bfsi().dataPtr(bnd_start),
			     ARLIM(bfsi().loVect()),ARLIM(bfsi().hiVect()),
			     bx.loVect(),bx.hiVect(),
			     mfi().dataPtr(mf_start),
			     ARLIM(mfi().loVect()),ARLIM(mfi().hiVect()),
			     &num_comp,hfine);
            }
        }
    }

    // now copy boundary values stored in ghost cells of fine
    // into bndry.  This only does something for physical boundaries,
    // we don't need to make it periodic aware
    for (OrientationIter fi; fi; ++fi) {
	bndry[fi()].copyFrom(mf,0,mf_start,bnd_start,num_comp);
    }
}


// (1) set bndry type and location of bndry value on each face of
//     each grid
// (2) set actual bndry value by:
//     (A) Interpolate from crse bndryRegister at crse/fine interface
//     (B) Copy from ghost region of MultiFab at physical bndry
//     (C) Copy from valid region of MultiFab at fine/fine interface
void
MCInterpBndryData::setBndryValues (const BndryRegister& crse,
				   int                  c_start,
				   const MultiFab&      fine,
				   int                  f_start,
				   int                  bnd_start,
				   int                  num_comp,
				   int                  ratio,
				   const Array<BCRec>&  bc)
{
    if (! bdfunc_set) bdfunc_init();

      // check that boxarrays are identical
    assert( grids.ready() );
    assert( grids == fine.boxArray() );

      // set bndry types and bclocs
    setBndryConds(bc, geom, ratio);

    // find the dx's
    const Real* h = geom.CellSize();

      // first interpolate from coarse to fine on bndry
    const BOX& fine_domain = geom.Domain();

      // mask turned off if covered by fine grid
    Real *derives = 0;
    int  tmplen = 0;
    for(ConstMultiFabIterator finemfi(fine); finemfi.isValid(); ++finemfi)
    {
        assert(grids[finemfi.index()] == finemfi.validbox());
        const Box &fine_bx = finemfi.validbox();
        BOX crse_bx = coarsen(fine_bx,ratio);
        const int* cblo = crse_bx.loVect();
        const int* cbhi = crse_bx.hiVect();
        int mxlen = crse_bx.longside() + 2;
        if (pow((double)mxlen,(double)BL_SPACEDIM-1) > tmplen) {
            delete derives;
            tmplen = mxlen;
#if (BL_SPACEDIM > 2)
	    tmplen *= mxlen;
#endif	    
            derives = new Real[tmplen*NUMDERIV];
        }
	const int* lo = fine_bx.loVect();
	const int* hi = fine_bx.hiVect();
	const FARRAYBOX& fine_grd = finemfi();
        const int* finelo = fine_grd.loVect();
        const int* finehi = fine_grd.hiVect();
        const Real* finedat = fine_grd.dataPtr(f_start);

	for (OrientationIter fi; fi; ++fi) {
	    Orientation face(fi());
	    int dir = face.coordDir();
	    // load up hfine with perpindicular h's
	    Real hfine[BL_SPACEDIM];
	    int kdir = 0;
	    for(int idir=0; idir<BL_SPACEDIM; ++idir){
	      if( idir == dir ) continue;
	      hfine[kdir++] = h[idir];
	    }
	    FARRAYBOX& bnd_fab = bndry[face][finemfi.index()];
	    const int* blo = bnd_fab.loVect();
	    const int* bhi = bnd_fab.hiVect();
	    Real* bdat = bnd_fab.dataPtr(bnd_start);

	    if (fine_bx[face] != fine_domain[face] ||
		geom.isPeriodic(dir)) {
		
		  // internal or periodic edge, interpolate from crse data
                const Mask& mask = masks[face][finemfi.index()];
                const int* mlo = mask.loVect();
                const int* mhi = mask.hiVect();
                const int* mdat = mask.dataPtr();

                const FARRAYBOX& crse_fab = crse[face][finemfi.index()];
                const int* clo = crse_fab.loVect();
                const int* chi = crse_fab.hiVect();
                const Real* cdat = crse_fab.dataPtr(c_start);

		int is_not_covered = BndryData::not_covered;
                bdfunc[face](bdat,ARLIM(blo),ARLIM(bhi),
			     lo,hi,ARLIM(cblo),ARLIM(cbhi),
                             &num_comp,&ratio,&is_not_covered,
			     mdat,ARLIM(mlo),ARLIM(mhi),
                             cdat,ARLIM(clo),ARLIM(chi),derives,hfine);
            } else {
		// this copies data from ghost region of grid, plus
		// computes derivative
		bdider[face](bdat,ARLIM(blo),ARLIM(bhi),
			     lo, hi,
			     finedat,ARLIM(finelo),ARLIM(finehi),
			     &num_comp,hfine);
	    }
	}
    }
    delete derives;

    // now copy boundary values stored in ghost cells of fine
    // into bndry.  This only does something for physical boundaries,
    // we don't need to make it periodic aware
    for (OrientationIter face; face; ++face) {
	bndry[face()].copyFrom(fine,0,f_start,bnd_start,num_comp);
    }
}

