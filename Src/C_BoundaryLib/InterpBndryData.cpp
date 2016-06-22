
#include <winstd.H>
#include <LO_BCTYPES.H>
#include <InterpBndryData.H>
#include <INTERPBNDRYDATA_F.H>

#if (BL_SPACEDIM == 1)
#define NUMDERIV 2
#endif

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

//
// For sliding parabolic interp in bdfuncs
//
int InterpBndryData::IBD_max_order_DEF = 3;

namespace
{
    bool initialized = false;

    BDInterpFunc* bdfunc[2*BL_SPACEDIM];
}

static
void
bdfunc_init ()
{
    const Orientation xloface(0,Orientation::low);
    const Orientation xhiface(0,Orientation::high);

    bdfunc[xloface] = FORT_BDINTERPXLO;
    bdfunc[xhiface] = FORT_BDINTERPXHI;

#if (BL_SPACEDIM > 1)
    const Orientation yloface(1,Orientation::low);
    const Orientation yhiface(1,Orientation::high);
    bdfunc[yloface] = FORT_BDINTERPYLO;
    bdfunc[yhiface] = FORT_BDINTERPYHI;
#endif

#if (BL_SPACEDIM > 2)
    const Orientation zloface(2,Orientation::low);
    const Orientation zhiface(2,Orientation::high);
    bdfunc[zloface] = FORT_BDINTERPZLO;
    bdfunc[zhiface] = FORT_BDINTERPZHI;
#endif
}

InterpBndryData::InterpBndryData ()
    :
    BndryData()
{}

InterpBndryData::InterpBndryData (const InterpBndryData& rhs)
    :
    BndryData(rhs)
{}

InterpBndryData&
InterpBndryData::operator= (const InterpBndryData& rhs)
{
    if (!(this == &rhs))
    {
        BndryData::operator=(rhs);
    }
    return *this;
}

InterpBndryData::InterpBndryData (const BoxArray& _grids,
                                  int             _ncomp,
                                  const Geometry& geom,
				  ParallelDescriptor::Color color)
    :
    BndryData(_grids,_ncomp,geom,color)
{}

InterpBndryData::~InterpBndryData () {}

void
InterpBndryData::setBndryConds (const BCRec& phys_bc,
                                int          ratio)
{

    const IntVect& ratio_vect = ratio * IntVect::TheUnitVector();
    setBndryConds(phys_bc, ratio_vect);
}

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
    BL_ASSERT(grids.size());
    BL_ASSERT(grids == mf.boxArray());

    IntVect ref_ratio = IntVect::TheUnitVector();

    for (int n = bnd_start; n < bnd_start+num_comp; ++n)
	setBndryConds(bc, ref_ratio, n);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        BL_ASSERT(grids[mfi.index()] == mfi.validbox());

        const Box& bx = grids[mfi.index()];

        for (OrientationIter fi; fi; ++fi)
        {
            const Orientation face = fi();

            if (bx[face]==geom.Domain()[face] && !geom.isPeriodic(face.coordDir()))
            {
                //
                // Physical bndry, copy from grid.
                //
                FArrayBox& bnd_fab = bndry[face][mfi];
                bnd_fab.copy(mf[mfi],mf_start,bnd_start,num_comp);
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
                                 const IntVect&  ratio,
                                 const BCRec&    bc,
                                 int             max_order)
{
    if (!initialized)
        bdfunc_init();
    //
    // Check that boxarrays are identical.
    //
    BL_ASSERT(grids.size());
    BL_ASSERT(grids == fine.boxArray());
    //
    // Set bndry types and bclocs.
    //
    for (int n = bnd_start; n < bnd_start+num_comp; ++n)
	setBndryConds(bc, ratio, n);
    //
    // First interpolate from coarse to fine on bndry.
    //
    const Box& fine_domain = geom.Domain();
    //
    // Mask turned off if covered by fine grid.
    //
    if (max_order==3 || max_order==1)
    {
        for (MFIter fine_mfi(fine); fine_mfi.isValid(); ++fine_mfi)
        {
            BL_ASSERT(grids[fine_mfi.index()] == fine_mfi.validbox());

            const Box&       fine_bx  = fine_mfi.validbox();
            const Box&       crse_bx  = BoxLib::coarsen(fine_bx,ratio);
            const int*       cblo     = crse_bx.loVect();
            const int*       cbhi     = crse_bx.hiVect();
            const int        mxlen    = crse_bx.longside() + 2;
            const int*       lo       = fine_bx.loVect();
            const int*       hi       = fine_bx.hiVect();
            const FArrayBox& fine_grd = fine[fine_mfi];
            const MaskTuple& msk      = masks[fine_mfi.index()];

#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < 2*BL_SPACEDIM; i++)
            {
                const int               dir  = ((i<BL_SPACEDIM) ? i : (i-BL_SPACEDIM));
                const Orientation::Side side = ((i<BL_SPACEDIM) ? Orientation::low : Orientation::high);

                const Orientation face(dir,side);

                if (fine_bx[face] != fine_domain[face] || geom.isPeriodic(dir))
                {
                    //
                    // Internal or periodic edge, interpolate from crse data.
                    //
                    Array<Real> derives(D_TERM(1,*mxlen,*mxlen)*NUMDERIV);

                    const Mask&      mask           = *(msk[face]);
                    const int*       mlo            = mask.loVect();
                    const int*       mhi            = mask.hiVect();
                    const int*       mdat           = mask.dataPtr();
                    const FArrayBox& crse_fab       = crse[face][fine_mfi];
                    const int*       clo            = crse_fab.loVect();
                    const int*       chi            = crse_fab.hiVect();
                    const Real*      cdat           = crse_fab.dataPtr(c_start);
                    FArrayBox&       bnd_fab        = bndry[face][fine_mfi];
                    const int*       blo            = bnd_fab.loVect();
                    const int*       bhi            = bnd_fab.hiVect();
                    Real*            bdat           = bnd_fab.dataPtr(bnd_start);
                    int              is_not_covered = BndryData::not_covered;
                    //
                    // The quadratic interp needs crse data in 2 grow cells tangential
                    // to face.  This checks to be sure the source data is large enough.
                    //
                    Box crsebnd = BoxLib::adjCell(crse_bx,face,1);

                    if (max_order == 3) 
                    {
                        for (int k=0;k<BL_SPACEDIM;k++)
                            if (k!=dir)
                                crsebnd.grow(k,2);
                        BL_ASSERT(crse_fab.box().contains(crsebnd));
                    }

                    bdfunc[face](bdat,ARLIM(blo),ARLIM(bhi),
                                 lo,hi,ARLIM(cblo),ARLIM(cbhi),
                                 &num_comp,ratio.getVect(),&is_not_covered,
                                 mdat,ARLIM(mlo),ARLIM(mhi),
                                 cdat,ARLIM(clo),ARLIM(chi),derives.dataPtr(),&max_order);
                }
                else
                {
                    //
                    // Physical bndry, copy from ghost region of corresponding grid
                    //
                    FArrayBox& bnd_fab = bndry[face][fine_mfi];
                    bnd_fab.copy(fine_grd,f_start,bnd_start,num_comp);
                }
            }
        }
    }
    else
    {
        BoxLib::Abort("InterpBndryData::setBndryValues supports only max_order=1 or 3");
    }
}

void
InterpBndryData::setBndryValues (::BndryRegister&  crse,
                                 int             c_start,
                                 const MultiFab& fine,
                                 int             f_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 int             ratio,
                                 const BCRec&    bc,
                                 int             max_order)
{
    const IntVect& ratio_vect = ratio * IntVect::TheUnitVector();
    setBndryValues(crse,c_start,fine,f_start,bnd_start,num_comp,ratio_vect,bc,max_order);
}
