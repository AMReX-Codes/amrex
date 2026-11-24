
#ifndef AMREX_INTERPBNDRYDATA_H_
#define AMREX_INTERPBNDRYDATA_H_
#include <AMReX_Config.H>

#include <AMReX_BCRec.H>
#include <AMReX_BndryData.H>
#include <AMReX_InterpBndryData_K.H>

/**
 * \file AMReX_InterpBndryData.H
 * \brief Boundary data specialization that interpolates coarse data into fine ghost cells.
 */

namespace amrex {

/// \cond DOXYGEN_IGNORE
namespace detail {
    // One boundary face of one box, for the fused GPU launch in setBndryValues.
    template <typename T>
    struct InterpBndryTag {
        Array4<T> bdry;
        Array4<T const> src;    // coarse data, or fine data on a physical boundary
        Array4<int const> mask;
        Box bx;
        int dir;                // direction of a coarse/fine face; -1 on a physical boundary

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Box const& box () const noexcept { return bx; }
    };
}
/// \endcond

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
     * \brief Sset bndry values at coarse/fine and physical boundaries
     *
     * \param crse      BndryRegister storing coarse level data
     * \param c_start   starting component of the coarse data
     * \param fine      MF containing physical boundary values
     * \param f_start   starting component of the fine data
     * \param bnd_start starting component in this InterpBndryData
     * \param num_comp  number of component
     * \param ratio     refinement ratio
     * \param max_order interpolation order
     * \param max_width max width of the interpolation stencil
     *
     * \note Order-3 interpolation reads the boundary mask at tangential offsets
     *       of `±ratio`, and the masks are allocated with a fixed tangential
     *       half-width (see `BndryDataT::NTangHalfWidth`).  Consequently
     *       `max_order == 3` requires `ratio <= 4` in every dimension; this is
     *       enforced at runtime.  `max_order == 1` has no such restriction.
     */
    void setBndryValues (BndryRegisterT<MF> const& crse, int c_start, const MF& fine, int f_start,
                         int bnd_start, int num_comp, const IntVect& ratio,
                         int max_order = IBD_max_order_DEF, int max_width = 2);

    /**
     * \brief Update boundary values at coarse/fine boundaries
     *
     * \param crse      BndryRegister storing coarse level data
     * \param c_start   starting component of the coarse data
     * \param bnd_start starting component in this InterpBndryData
     * \param num_comp  number of component
     * \param ratio     refinement ratio
     * \param max_order interpolation order
     * \param max_width max width of the interpolation stencil
     */
    void updateBndryValues (BndryRegisterT<MF>& crse, int c_start, int bnd_start, int num_comp,
                            const IntVect& ratio, int max_order = IBD_max_order_DEF, int max_width = 2);

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
                                      int max_order, int max_width)
{
    AMREX_ASSERT(this->grids == fine.boxArray());

    const Box& fine_domain = this->geom.Domain();

    // Order-3 interpolation reads the boundary mask at tangential offsets of
    // ±ratio.  The masks are built with a fixed tangential half-width
    // (BndryDataT::NTangHalfWidth == 5, i.e. ref_ratio + 1), so ratios above 4
    // would walk past the allocated mask halo.  Order-1 does not touch the mask.
    if (max_order == 3) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ratio.allLE(4),
            "InterpBndryDataT<MF>::setBndryValues: max_order=3 requires refinement ratio <= 4");
    }

#ifdef AMREX_USE_GPU
    if ((max_order==3 || max_order==1) && Gpu::inLaunchRegion())
    {
        // One launch for all faces of all boxes.
        using Tag = detail::InterpBndryTag<value_type>;
        Vector<Tag> tags;
        tags.reserve(std::size_t(fine.local_size())*std::size_t(2*AMREX_SPACEDIM));
        for (MFIter mfi(fine,MFItInfo().DisableDeviceSync()); mfi.isValid(); ++mfi)
        {
            const Box& fine_bx  = mfi.validbox();
            for (OrientationIter fi; fi; ++fi) {
                const Orientation face = fi();
                const int dir = face.coordDir();
                if (fine_bx[face] != fine_domain[face] || this->geom.isPeriodic(dir))
                {
                    tags.push_back(Tag{.bdry = this->bndry[face].array(mfi),
                                       .src  = crse[face].const_array(mfi),
                                       .mask = (max_order == 3) ? this->masks[face].const_array(mfi)
                                                                : Array4<int const>{},
                                       .bx   = this->bndry[face][mfi].box(),
                                       .dir  = dir});
                }
                else if (fine.defined(mfi)) {
                    auto      & bnd_fab = this->bndry[face][mfi];
                    auto const& src_fab = fine[mfi];
                    tags.push_back(Tag{.bdry = bnd_fab.array(),
                                       .src  = src_fab.const_array(),
                                       .mask = Array4<int const>{},
                                       .bx   = bnd_fab.box() & src_fab.box(),
                                       .dir  = -1});
                }
            }
        }

        const auto rr = ratio.dim3();
        const int is_not_covered = BndryData::not_covered;
        ParallelFor(tags, num_comp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n, Tag const& tag) noexcept
        {
            if (tag.dir < 0) {
                tag.bdry(i,j,k,n+bnd_start) = tag.src(i,j,k,n+f_start);
            } else if (max_order == 1) {
                interpbndrydata_o1(i,j,k,n,tag.bdry,bnd_start,tag.src,c_start,rr);
            } else if (tag.dir == 0) {
                interpbndrydata_x_o3(i,j,k,n,tag.bdry,bnd_start,tag.src,c_start,rr,
                                     tag.mask, is_not_covered, max_width);
            }
#if (AMREX_SPACEDIM >= 2)
            else if (tag.dir == 1) {
                interpbndrydata_y_o3(i,j,k,n,tag.bdry,bnd_start,tag.src,c_start,rr,
                                     tag.mask, is_not_covered, max_width);
            }
#if (AMREX_SPACEDIM == 3)
            else {
                interpbndrydata_z_o3(i,j,k,n,tag.bdry,bnd_start,tag.src,c_start,rr,
                                     tag.mask, is_not_covered, max_width);
            }
#endif
#endif
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    }
    else
#endif
    if (max_order==3 || max_order==1)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(fine,MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
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
                        amrex::LoopConcurrentOnCpu(b, num_comp, [&] (int it, int jt, int kt, int n)
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
                            amrex::LoopConcurrentOnCpu(b, num_comp, [&] (int it, int jt, int kt, int n)
                            {
                                interpbndrydata_x_o3(it,jt,kt,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr,
                                                     mask_array, is_not_covered, max_width);
                            });
                            break;
                        }
#if (AMREX_SPACEDIM >= 2)
                        case 1:
                        {
                            amrex::LoopConcurrentOnCpu(b, num_comp, [&] (int it, int jt, int kt, int n)
                            {
                                interpbndrydata_y_o3(it,jt,kt,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr,
                                                     mask_array, is_not_covered, max_width);
                            });
                            break;
                        }
#if (AMREX_SPACEDIM == 3)
                        case 2:
                        {
                            amrex::LoopConcurrentOnCpu(b, num_comp, [&] (int it, int jt, int kt, int n)
                            {
                                interpbndrydata_z_o3(it,jt,kt,n,bdry_array,bnd_start,
                                                     crse_array,c_start,rr,
                                                     mask_array, is_not_covered, max_width);
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
                    amrex::LoopConcurrentOnCpu(b, num_comp, [&] (int ii, int jj, int kk, int nn)
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
                                         int num_comp, const IntVect& ratio, int max_order, int max_width)
{
    MF foo(this->grids, this->bndry[0].DistributionMap(), 1, num_comp, MFInfo().SetAlloc(false));
    setBndryValues(crse, c_start, foo, 0, bnd_start, num_comp, ratio, max_order, max_width);
}

template <typename MF>
void
InterpBndryDataT<MF>::setHomogValues ()
{
    this->setVal(0);
}

using InterpBndryData = InterpBndryDataT<MultiFab>;
using fInterpBndryData = InterpBndryDataT<fMultiFab>;

}

#endif /*AMREX_INTERPBNDRYDATA_H_*/
