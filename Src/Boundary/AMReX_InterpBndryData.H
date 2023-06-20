
#ifndef AMREX_INTERPBNDRYDATA_H_
#define AMREX_INTERPBNDRYDATA_H_
#include <AMReX_Config.H>

#include <AMReX_BCRec.H>
#include <AMReX_BndryData.H>
#include <AMReX_InterpBndryData_K.H>

namespace amrex {

/**
 * \brief An InterpBndryData object adds to a BndryData object the ability to
        manipulate and set the data stored in the boundary cells.

        The InterpBndryData class is a class derived from
        BndryData.  It is intended to provide a more physical method for
        filling boundary-related data.  Boundary values in a BndryData object
        are stored in FabSets around each grid in the domain, and the
        InterpBndryData class provides a mechanism for filling these FabSets,
        consistent with AMR-like numerical discretizations.  When asked to
        set its boundary values, an InterpBndryData object:

        Fills with physical boundary values if the FAB is on the
        domain boundary - the corresponding values are presumed to be
        stored in the ghost cells of a MultiFab given to the boundary filling
        routine

        Fills on intersection with data from the VALID region of the
        input MultiFab, ensuring that adjacent FABs use consistent data at
        their intersection, and otherwise,

        Fills with values interpolated from a coarser FAB that
        bounds the cells that do not meet the above two criteria
*/
template <typename MF>
class InterpBndryDataT
    :
    public BndryDataT<MF>
{
public:

    using value_type = typename MF::value_type;

    /**
    * \brief default constructor
    */
    InterpBndryDataT () noexcept = default;

    /**
    * \brief constructor for given BoxArray, etc
    *
    * \param _grids
    * \param _dmap
    * \param _ncomp
    * \param geom
    */
    InterpBndryDataT (const BoxArray& _grids,
                      const DistributionMapping& _dmap,
                      int             _ncomp,
                      const Geometry& geom);

    ~InterpBndryDataT () = default;

    InterpBndryDataT (const InterpBndryDataT<MF>& rhs) = delete;
    InterpBndryDataT (InterpBndryDataT<MF>&& rhs) = delete;
    InterpBndryDataT<MF>& operator= (const InterpBndryDataT<MF>& rhs) = delete;
    InterpBndryDataT<MF>& operator= (InterpBndryDataT<MF>&& rhs) = delete;

    /**
     * \brief Set bndry values at physical boundaries
     *
     * \param mf        MF containing physical boundary values
     * \param mf_start  starting component of the data in MF
     * \param bnd_start starting component in this boundary register
     * \param num_comp  number of components
     */
    void setPhysBndryValues (const MF& mf, int mf_start, int bnd_start, int num_comp);

    /**
     * \briefSset bndry values at coarse/fine and physical boundaries
     *
     * \param crse      BndryRegister storing coarse level data
     * \param c_start   starting component of the coarse data
     * \param fine      MF containing physical bounadry values
     * \param f_start   starting component of the fine data
     * \param bnd_start starting component in this InterpBndryData
     * \param num_comp  number of component
     * \param ratio     refinement ratio
     * \param max_order interpolation order
     */
    void setBndryValues (BndryRegisterT<MF> const& crse, int c_start, const MF& fine, int f_start,
                         int bnd_start, int num_comp, const IntVect& ratio,
                         int max_order = IBD_max_order_DEF);

    /**
     * \brief Update boundary values at coarse/fine boundaries
     *
     * \param crse      BndryRegister storing coarse level data
     * \param c_start   starting component of the coarse data
     * \param bnd_start starting component in this InterpBndryData
     * \param num_comp  number of component
     * \param ratio     refinement ratio
     * \param max_order interpolation order
     */
    void updateBndryValues (BndryRegisterT<MF>& crse, int c_start, int bnd_start, int num_comp,
                            const IntVect& ratio, int max_order = IBD_max_order_DEF);

    //! Set boundary values to zero
    void setHomogValues ();

    // For sliding parabolic interp in bdfuncs
    static constexpr int IBD_max_order_DEF = 3;

};

template <typename MF>
InterpBndryDataT<MF>::InterpBndryDataT (const BoxArray& _grids,
                                        const DistributionMapping& _dmap,
                                        int             _ncomp,
                                        const Geometry& _geom)
    :
    BndryDataT<MF>(_grids,_dmap,_ncomp,_geom)
{}

template <typename MF>
void
InterpBndryDataT<MF>::setPhysBndryValues (const MF& mf, int mf_start, int bnd_start, int num_comp)
{
    AMREX_ASSERT(this->grids == mf.boxArray());

    const Box& fine_domain = this->geom.Domain();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf,MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        for (OrientationIter fi; fi; ++fi) {
            const Orientation face = fi();
            if (bx[face] == fine_domain[face] && !(this->geom.isPeriodic(face.coordDir())))
            {
                // Physical bndry, copy from grid.
                auto      & bnd_fab = this->bndry[face][mfi];
                auto const& src_fab = mf[mfi];
                auto const& bnd_array = bnd_fab.array();
                auto const& src_array = src_fab.const_array();
                const Box& b = src_fab.box() & bnd_fab.box();
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( b, num_comp, i, j, k, n,
                {
                    bnd_array(i,j,k,n+bnd_start) = src_array(i,j,k,n+mf_start);
                });
            }
        }
    }
}

template <typename MF>
void
InterpBndryDataT<MF>::setBndryValues (BndryRegisterT<MF> const& crse, int c_start, const MF& fine,
                                      int f_start, int bnd_start, int num_comp, const IntVect& ratio,
                                      int max_order)
{
    AMREX_ASSERT(this->grids == fine.boxArray());

    const Box& fine_domain = this->geom.Domain();

    if (max_order==3 || max_order==1)
    {
        MFItInfo info;
        if (Gpu::notInLaunchRegion()) info.SetDynamic(true);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(fine,info); mfi.isValid(); ++mfi)
        {
            const Box& fine_bx  = mfi.validbox();
            for (OrientationIter fi; fi; ++fi) {
                const Orientation face = fi();
                const int dir = face.coordDir();
                if (fine_bx[face] != fine_domain[face] || this->geom.isPeriodic(dir))
                { // coarse/fine boundary: interpolate from coarse data stored in BndryRegister crse
                    auto const& crse_array = crse[face].const_array(mfi);
                    auto const& bdry_array = this->bndry[face].array(mfi);
                    Box const& b = this->bndry[face][mfi].box();
                    const auto rr = ratio.dim3();

                    if (max_order == 1)
                    {
                        AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, it, jt, kt, n,
                        {
                            interpbndrydata_o1(it,jt,kt,n,bdry_array,bnd_start,
                                               crse_array,c_start,rr);
                        });
                    }
                    else
                    {
                        auto const& mask_array = this->masks[face].const_array(mfi);
                        int is_not_covered = BndryData::not_covered;
                        switch (dir) {
                        case 0:
                        {
                            AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, it, jt, kt, n,
                            {
                                interpbndrydata_x_o3(it,jt,kt,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr,
                                                     mask_array, is_not_covered);
                            });
                            break;
                        }
#if (AMREX_SPACEDIM >= 2)
                        case 1:
                        {
                            AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, it, jt, kt, n,
                            {
                                interpbndrydata_y_o3(it,jt,kt,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr,
                                                     mask_array, is_not_covered);
                            });
                            break;
                        }
#if (AMREX_SPACEDIM == 3)
                        case 2:
                        {
                            AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, it, jt, kt, n,
                            {
                                interpbndrydata_z_o3(it,jt,kt,n,bdry_array,bnd_start,
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
                else if (fine.defined(mfi)) {
                    // Physical bndry
                    auto      & bnd_fab = this->bndry[face][mfi];
                    auto const& src_fab = fine[mfi];
                    auto const& bnd_array = bnd_fab.array();
                    auto const& src_array = src_fab.const_array();
                    const Box& b = bnd_fab.box() & src_fab.box();
                    AMREX_HOST_DEVICE_FOR_4D ( b, num_comp, ii, jj, kk, nn,
                    {
                        bnd_array(ii,jj,kk,nn+bnd_start) = src_array(ii,jj,kk,nn+f_start);
                    });
                }
            }
        }
    }
    else
    {
        amrex::Abort("InterpBndryDataT<MF>::setBndryValues supports only max_order=1 or 3");
    }
}

template <typename MF>
void
InterpBndryDataT<MF>::updateBndryValues (BndryRegisterT<MF>& crse, int c_start, int bnd_start,
                                         int num_comp, const IntVect& ratio, int max_order)
{
    MF foo(this->grids, this->bndry[0].DistributionMap(), 1, num_comp, MFInfo().SetAlloc(false));
    setBndryValues(crse, c_start, foo, 0, bnd_start, num_comp, ratio, max_order);
}

template <typename MF>
void
InterpBndryDataT<MF>::setHomogValues ()
{
    this->setVal(0.);
}

using InterpBndryData = InterpBndryDataT<MultiFab>;
using fInterpBndryData = InterpBndryDataT<fMultiFab>;

}

#endif /*AMREX_INTERPBNDRYDATA_H_*/
