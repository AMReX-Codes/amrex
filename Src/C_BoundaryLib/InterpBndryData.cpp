//BL_COPYRIGHT_NOTICE

//
// $Id: InterpBndryData.cpp,v 1.10 1999-06-28 20:37:04 car Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cmath>
#else
#include <math.h>
#endif

#include <LO_BCTYPES.H>
#include <InterpBndryData.H>
#include <INTERPBNDRYDATA_F.H>

static BDInterpFunc* bdfunc[2*BL_SPACEDIM];

static int bdfunc_set = 0;

static
void
bdfunc_init ()
{
    Orientation xloface(0,Orientation::low);
    Orientation xhiface(0,Orientation::high);

    bdfunc[xloface] = FORT_BDINTERPXLO;
    bdfunc[xhiface] = FORT_BDINTERPXHI;

#if (BL_SPACEDIM > 1)
    Orientation yloface(1,Orientation::low);
    Orientation yhiface(1,Orientation::high);
    bdfunc[yloface] = FORT_BDINTERPYLO;
    bdfunc[yhiface] = FORT_BDINTERPYHI;
#endif

#if (BL_SPACEDIM > 2)
    Orientation zloface(2,Orientation::low);
    Orientation zhiface(2,Orientation::high);
    bdfunc[zloface] = FORT_BDINTERPZLO;
    bdfunc[zhiface] = FORT_BDINTERPZHI;
#endif
}

#if (BL_SPACEDIM == 2)
#define NUMDERIV 2
#endif

#if (BL_SPACEDIM == 3)
#define NUMDERIV 5
#endif

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();
#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

InterpBndryData::InterpBndryData (const BoxArray& _grids,
                                  int             _ncomp,
                                  const Geometry& geom)
    :
    BndryData(_grids,_ncomp,geom)
{}

//
// At the coarsest level the bndry values are taken from adjacent grids.
//

void
InterpBndryData::setBndryValues (const MultiFab& mf,
                                 int             mf_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 const BCRec&    bc)
{
    //
    // Check that boxarrays are identical.
    //
    BL_ASSERT(grids.ready());
    BL_ASSERT(grids == mf.boxArray());

    IntVect ref_ratio = IntVect::TheUnitVector();
    for (int n=bnd_start; n<bnd_start+num_comp; ++n)
	setBndryConds(bc, ref_ratio, n);

    for (ConstMultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
    {
        BL_ASSERT(grids[mfi.index()] == mfi.validbox());

        const Box& bx = grids[mfi.index()];

        for (OrientationIter fi; fi; ++fi)
        {
            Orientation face(fi());

            if (bx[face]==geom.Domain()[face] && !geom.isPeriodic(face.coordDir()))
            {
                //
                // Physical bndry, copy from grid.
                //
                FArrayBox& bnd_fab = bndry[face][mfi.index()];
                bnd_fab.copy(mf[mfi.index()],mf_start,bnd_start,num_comp);
            }
        }
    }
}

//
// (1) set bndry type and location of bndry value on each face of
//     each grid
// (2) set actual bndry value by:
//     (A) Interpolate from crse bndryRegister at crse/fine interface
//     (B) Copy from ghost region of MultiFab at physical bndry
//

void
InterpBndryData::setBndryValues (::BndryRegister& crse,
                                 int             c_start,
                                 const MultiFab& fine,
                                 int             f_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 IntVect&        ratio,
                                 const BCRec&    bc)
{
    if (!bdfunc_set)
        bdfunc_init();
    //
    // Check that boxarrays are identical.
    //
    BL_ASSERT(grids.ready());
    BL_ASSERT(grids == fine.boxArray());
    //
    // Set bndry types and bclocs.
    //
    setBndryConds(bc, ratio);
    //
    // First interpolate from coarse to fine on bndry.
    //
    const Box& fine_domain = geom.Domain();
    //
    // Mask turned off if covered by fine grid.
    //
    Real* derives = 0;
    int tmplen    = 0;

    for (ConstMultiFabIterator fine_mfi(fine); fine_mfi.isValid(); ++fine_mfi)
    {
        BL_ASSERT(grids[fine_mfi.index()] == fine_mfi.validbox());

        const Box& fine_bx = fine_mfi.validbox();
        Box crse_bx        = coarsen(fine_bx,ratio);
        const int* cblo    = crse_bx.loVect();
        const int* cbhi    = crse_bx.hiVect();
        int mxlen          = crse_bx.longside() + 2;

        if (pow(float(mxlen), float(BL_SPACEDIM-1)) > tmplen)
        {
            delete derives;
            tmplen = mxlen;
#if (BL_SPACEDIM > 2)
            tmplen *= mxlen;
#endif            
            derives = new Real[tmplen*NUMDERIV];
        }
        const int* lo             = fine_bx.loVect();
        const int* hi             = fine_bx.hiVect();
        const FArrayBox &fine_grd = fine_mfi();
        const int *finelo         = fine_grd.loVect();
        const int *finehi         = fine_grd.hiVect();
        const Real *finedat       = fine_grd.dataPtr(f_start);

        for (OrientationIter fi; fi; ++fi)
        {
            Orientation face(fi());
            int dir = face.coordDir();
            if (fine_bx[face] != fine_domain[face] || geom.isPeriodic(dir))
            {
                //
                // Internal or periodic edge, interpolate from crse data.
                //
                const Mask& mask          = masks[face][fine_mfi.index()];
                const int* mlo            = mask.loVect();
                const int* mhi            = mask.hiVect();
                const int* mdat           = mask.dataPtr();
                const FArrayBox& crse_fab = crse[face][fine_mfi.index()];
                const int* clo            = crse_fab.loVect();
                const int* chi            = crse_fab.hiVect();
                const Real* cdat          = crse_fab.dataPtr(c_start);
                FArrayBox& bnd_fab        = bndry[face][fine_mfi.index()];
                const int* blo            = bnd_fab.loVect();
                const int* bhi            = bnd_fab.hiVect();
                Real* bdat                = bnd_fab.dataPtr(bnd_start);
                int is_not_covered        = BndryData::not_covered;

                bdfunc[face](bdat,ARLIM(blo),ARLIM(bhi),
                             lo,hi,ARLIM(cblo),ARLIM(cbhi),
                             &num_comp,ratio.getVect(),&is_not_covered,
                             mdat,ARLIM(mlo),ARLIM(mhi),
                             cdat,ARLIM(clo),ARLIM(chi),derives);
            }
            else
            {
                //
                // Physical bndry, copy from ghost region of corresponding grid
                //
                FArrayBox &bnd_fab = bndry[face][fine_mfi.index()];
                bnd_fab.copy(fine_grd,f_start,bnd_start,num_comp);
            }
        }
    }
    delete derives;
}
