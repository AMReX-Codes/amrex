
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_InterpBndryData.H>
#include <AMReX_INTERPBNDRYDATA_F.H>

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

namespace amrex {

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
				  const DistributionMapping& _dmap,
                                  int             _ncomp,
                                  const Geometry& _geom)
    :
    BndryData(_grids,_dmap,_ncomp,_geom)
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
    setBndryValues(mf, mf_start, bnd_start, num_comp, IntVect::TheUnitVector(), bc);
}

void
InterpBndryData::setBndryValues (const MultiFab& mf,
                                 int             mf_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 int             ref_ratio,
                                 const BCRec&    bc)
{
    setBndryValues(mf, mf_start, bnd_start, num_comp, IntVect{ref_ratio}, bc);
}

void
InterpBndryData::setBndryValues (const MultiFab& mf,
                                 int             mf_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 const IntVect&  ref_ratio,
                                 const BCRec&    bc)
{
    //
    // Check that boxarrays are identical.
    //
    BL_ASSERT(grids.size());
    BL_ASSERT(grids == mf.boxArray());

    for (int n = bnd_start; n < bnd_start+num_comp; ++n) {
	setBndryConds(bc, ref_ratio, n);
    }

    // TODO: tiling - wqz
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf,MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

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
InterpBndryData::setBndryValues (BndryRegister& crse,
                                 int             c_start,
                                 const MultiFab& fine,
                                 int             f_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 const IntVect&  ratio,
                                 const BCRec&    bc,
                                 int             max_order)
{
    if (!initialized) {
        bdfunc_init();
    }

    BndryValuesDoIt (crse, c_start, &fine, f_start, bnd_start, num_comp, ratio, &bc, max_order);
}

void
InterpBndryData::BndryValuesDoIt (BndryRegister&  crse,
                                  int             c_start,
                                  const MultiFab* fine,
                                  int             f_start,
                                  int             bnd_start,
                                  int             num_comp,
                                  const IntVect&  ratio, 
                                  const BCRec*    bc,
                                  int             max_order)
{
    //
    // Check that boxarrays are identical.
    //
    BL_ASSERT(grids.size());
    BL_ASSERT(fine == nullptr || grids == fine->boxArray());
    //
    // Set bndry types and bclocs.
    //
    if (bc != nullptr) {
        for (int n = bnd_start; n < bnd_start+num_comp; ++n) {
            setBndryConds(*bc, ratio, n);
        }
    }
    //
    // First interpolate from coarse to fine on bndry.
    //
    const Box& fine_domain = geom.Domain();
    //
    // Mask turned off if covered by fine grid.
    //
    if (max_order==3 || max_order==1)
    {
        MultiFab foo(grids,bndry[0].DistributionMap(), 1, 0, MFInfo().SetAlloc(false), FArrayBoxFactory());

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
        Vector<Real> derives;

        for (MFIter mfi(foo,MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
        {
            BL_ASSERT(grids[mfi.index()] == mfi.validbox());

            const Box&       fine_bx  = mfi.validbox();
            const Box&       crse_bx  = amrex::coarsen(fine_bx,ratio);
            const int*       cblo     = crse_bx.loVect();
            const int*       cbhi     = crse_bx.hiVect();
            const int        mxlen    = crse_bx.longside() + 2;
            const int*       lo       = fine_bx.loVect();
            const int*       hi       = fine_bx.hiVect();

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
                    derives.resize(AMREX_D_TERM(1,*mxlen,*mxlen)*NUMDERIV);

                    const Mask&      mask           = masks[face][mfi];
                    const int*       mlo            = mask.loVect();
                    const int*       mhi            = mask.hiVect();
                    const int*       mdat           = mask.dataPtr();
                    const FArrayBox& crse_fab       = crse[face][mfi];
                    const int*       clo            = crse_fab.loVect();
                    const int*       chi            = crse_fab.hiVect();
                    const Real*      cdat           = crse_fab.dataPtr(c_start);
                    FArrayBox&       bnd_fab        = bndry[face][mfi];
                    const int*       blo            = bnd_fab.loVect();
                    const int*       bhi            = bnd_fab.hiVect();
                    Real*            bdat           = bnd_fab.dataPtr(bnd_start);
                    int              is_not_covered = BndryData::not_covered;
                    //
                    // The quadratic interp needs crse data in 2 grow cells tangential
                    // to face.  This checks to be sure the source data is large enough.
                    //
                    Box crsebnd = amrex::adjCell(crse_bx,face,1);

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
                else if (fine != nullptr)
                {
                    //
                    // Physical bndry, copy from ghost region of corresponding grid
                    //
                    const FArrayBox& fine_grd = (*fine)[mfi];
                    FArrayBox& bnd_fab = bndry[face][mfi];
                    bnd_fab.copy(fine_grd,f_start,bnd_start,num_comp);
                }
            }
        }
        }
    }
    else
    {
        amrex::Abort("InterpBndryData::setBndryValues supports only max_order=1 or 3");
    }
}

void
InterpBndryData::setBndryValues (BndryRegister&  crse,
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

void
InterpBndryData::updateBndryValues (BndryRegister& crse, int c_start, int bnd_start, int num_comp,
                                    const IntVect& ratio, int max_order)
{
    BndryValuesDoIt(crse, c_start, nullptr, 0, bnd_start, num_comp, ratio, nullptr, max_order);
}

void
InterpBndryData::updateBndryValues (BndryRegister& crse, int c_start, int bnd_start, int num_comp,
                                    int ratio, int max_order)
{
    updateBndryValues(crse, c_start, bnd_start, num_comp, IntVect{ratio}, max_order);
}

void
InterpBndryData::setHomogValues ()
{
    for (OrientationIter fi; fi; ++fi)
    {
        const Orientation face  = fi();
        
        for (FabSetIter fsi(bndry[face]); fsi.isValid(); ++fsi)
        {
            bndry[face][fsi].setVal(0.);
        }
    }
}

}
