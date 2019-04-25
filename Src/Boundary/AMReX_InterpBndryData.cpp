
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_InterpBndryData.H>
#include <AMReX_InterpBndryData_K.H>

namespace amrex {

//
// For sliding parabolic interp in bdfuncs
//
int InterpBndryData::IBD_max_order_DEF = 3;

InterpBndryData::InterpBndryData () noexcept
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

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
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
                auto bnd_fab = bndry[face].fabHostPtr(mfi);
                auto src_fab = mf.fabHostPtr(mfi);
                auto bnd_array = bnd_fab->array();
                auto const src_array = src_fab->array();
                const Box& b = src_fab->box() & bnd_fab->box();
                AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, i, j, k, n,
                {
                    bnd_array(i,j,k,n+bnd_start) = src_array(i,j,k,n+mf_start);
                });
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
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
        Vector<Real> derives;

        for (MFIter mfi(foo,MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
        {
            BL_ASSERT(grids[mfi.index()] == mfi.validbox());

            const Box&       fine_bx  = mfi.validbox();
            const Box&       crse_bx  = amrex::coarsen(fine_bx,ratio);

            for (int i = 0; i < 2*AMREX_SPACEDIM; i++)
            {
                const int               dir  = ((i<AMREX_SPACEDIM) ? i : (i-AMREX_SPACEDIM));
                const Orientation::Side side = ((i<AMREX_SPACEDIM) ? Orientation::low : Orientation::high);

                const Orientation face(dir,side);

                if (fine_bx[face] != fine_domain[face] || geom.isPeriodic(dir))
                {
                    auto const crse_array = crse[face].array(mfi);
                    auto       bdry_array = bndry[face].array(mfi);
                    const int dirside = dir*2 + side;
                    const auto rr = ratio.dim3();
                    //
                    // Internal or periodic edge, interpolate from crse data.
                    //
                    if (max_order == 1)
                    {
                        switch (dirside) {
                        case 0:
                        {
                            const Box& b = amrex::adjCellLo(crse_bx, 0);
                            AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, ic, jc, kc, n,
                            {
                                interpbndrydata_x_o1(1,ic,jc,kc,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr);
                            });
                            break;
                        }
                        case 1:
                        {
                            const Box& b = amrex::adjCellHi(crse_bx, 0);
                            AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, ic, jc, kc, n,
                            {
                                interpbndrydata_x_o1(0,ic,jc,kc,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr);
                            });
                            break;
                        }
#if (AMREX_SPACEDIM >= 2)
                        case 2:
                        {
                            const Box& b = amrex::adjCellLo(crse_bx, 1);
                            AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, ic, jc, kc, n,
                            {
                                interpbndrydata_y_o1(1,ic,jc,kc,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr);
                            });
                            break;
                        }
                        case 3:
                        {
                            const Box& b = amrex::adjCellHi(crse_bx, 1);
                            AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, ic, jc, kc, n,
                            {
                                interpbndrydata_y_o1(0,ic,jc,kc,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr);
                            });
                            break;
                        }
#if (AMREX_SPACEDIM == 3)
                        case 4:
                        {
                            const Box& b = amrex::adjCellLo(crse_bx, 2);
                            AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, ic, jc, kc, n,
                            {
                                interpbndrydata_z_o1(1,ic,jc,kc,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr);
                            });
                            break;
                        }
                        case 5:
                        {
                            const Box& b = amrex::adjCellHi(crse_bx, 2);
                            AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, ic, jc, kc, n,
                            {
                                interpbndrydata_z_o1(0,ic,jc,kc,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr);
                            });
                            break;
                        }
#endif
#endif
                        default: {}
                        }
                    }
                    else
                    {
                        auto const mask_array = masks[face].array(mfi);
                        int is_not_covered = BndryData::not_covered;
                        switch (dirside) {
                        case 0:
                        {
                            const Box& b = amrex::adjCellLo(crse_bx, 0);
                            AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, ic, jc, kc, n,
                            {
                                interpbndrydata_x_o3(1,ic,jc,kc,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr,
                                                     mask_array, is_not_covered);
                            });
                            break;
                        }
                        case 1:
                        {
                            const Box& b = amrex::adjCellHi(crse_bx, 0);
                            AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, ic, jc, kc, n,
                            {
                                interpbndrydata_x_o3(0,ic,jc,kc,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr,
                                                     mask_array, is_not_covered);
                            });
                            break;
                        }
#if (AMREX_SPACEDIM >= 2)
                        case 2:
                        {
                            const Box& b = amrex::adjCellLo(crse_bx, 1);
                            AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, ic, jc, kc, n,
                            {
                                interpbndrydata_y_o3(1,ic,jc,kc,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr,
                                                     mask_array, is_not_covered);
                            });
                            break;
                        }
                        case 3:
                        {
                            const Box& b = amrex::adjCellHi(crse_bx, 1);
                            AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, ic, jc, kc, n,
                            {
                                interpbndrydata_y_o3(0,ic,jc,kc,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr,
                                                     mask_array, is_not_covered);
                            });
                            break;
                        }
#if (AMREX_SPACEDIM == 3)
                        case 4:
                        {
                            const Box& b = amrex::adjCellLo(crse_bx, 2);
                            AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, ic, jc, kc, n,
                            {
                                interpbndrydata_z_o3(1,ic,jc,kc,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr,
                                                     mask_array, is_not_covered);
                            });
                            break;
                        }
                        case 5:
                        {
                            const Box& b = amrex::adjCellHi(crse_bx, 2);
                            AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, ic, jc, kc, n,
                            {
                                interpbndrydata_z_o3(0,ic,jc,kc,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr,
                                                     mask_array, is_not_covered);
                            });
                            break;
                        }
#endif
#endif
                        default: {}
                        }
                    }
                }
                else if (fine != nullptr)
                {
                    //
                    // Physical bndry, copy from ghost region of corresponding grid
                    //
                    auto bnd_fab = bndry[face].fabHostPtr(mfi);
                    auto src_fab = fine->fabHostPtr(mfi);
                    auto bnd_array = bnd_fab->array();
                    auto const src_array = src_fab->array();
                    const Box& b = bnd_fab->box() & src_fab->box();
                    AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, ii, jj, kk, nn,
                    {
                        bnd_array(ii,jj,kk,nn+bnd_start) = src_array(ii,jj,kk,nn+f_start);
                    });
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
    setVal(0.);
}

}
