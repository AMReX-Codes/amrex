#ifndef AMREX_MULTIMASK_H_
#define AMREX_MULTIMASK_H_
#include <AMReX_Config.H>

#include <AMReX_Mask.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>

/**
 * \file AMReX_MultiMask.H
 * \brief MultiFab-style container that stores Mask objects per grid.
 */

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

namespace amrex {

class MultiMaskIter;

class MultiMask
{
public:
    //! Construct an empty MultiMask; call define() before use.
    MultiMask () { ; }

    //! Construct masks over \p ba with mapping \p dm and \p ncomp components.
    MultiMask (const BoxArray& ba, const DistributionMapping& dm, int ncomp);

    /**
     * \brief Construct masks surrounding coarse/fine interfaces.
     *
     * \param regba  Regular BoxArray to cover.
     * \param dm     Distribution mapping.
     * \param geom   Geometry describing the physical domain.
     * \param face   Orientation for which to build masks.
     * \param in_rad Number of inward cells.
     * \param out_rad Number of outward cells.
     * \param extent_rad Tangential extent.
     * \param ncomp  Number of components.
     * \param initval Whether to initialize mask values immediately.
     */
    MultiMask (const BoxArray& regba, const DistributionMapping& dm, const Geometry& geom,
               Orientation face, int in_rad, int out_rad, int extent_rad, int ncomp, bool initval);

    ~MultiMask () = default;

    MultiMask (MultiMask&& rhs) noexcept = default;

    MultiMask (const MultiMask& rhs) = delete;
    MultiMask& operator= (const MultiMask& rhs) = delete;
    MultiMask& operator= (MultiMask&& rhs) = delete;

    //! Define masks over \p ba with mapping \p dm and \p ncomp components.
    void define (const BoxArray& ba, const DistributionMapping& dm, int ncomp);

    //! Define masks for a coarse/fine interface region described by \p regba, \p dm, \p geom, \p face, and stencil widths \p in_rad/\p out_rad/\p extent_rad.
    void define (const BoxArray& regba, const DistributionMapping& dm, const Geometry& geom,
                 Orientation face, int in_rad, int out_rad, int extent_rad, int ncomp, bool initval);

    //! Mutable access to the mask referenced by iterator \p mfi.
    Mask& operator[] (const MFIter& mfi) noexcept { return m_fa[mfi]; }

    //! Const access to the mask referenced by iterator \p mfi.
    const Mask& operator[] (const MFIter& mfi) const  noexcept { return m_fa[mfi]; }

    //! Return an Array4 view (const) for iterator \p mfi.
    Array4<int const> array (const MFIter& mfi) const noexcept { return m_fa.array(mfi); }

    //! Return an Array4 view (mutable) for iterator \p mfi.
    Array4<int      > array (const MFIter& mfi)       noexcept { return m_fa.array(mfi); }

    //! Return a const Array4 view for iterator \p mfi.
    Array4<int const> const_array (const MFIter& mfi) const noexcept {
        return m_fa.const_array(mfi);
    }

    //! Return const multi-array views of every mask.
    MultiArray4<int const> arrays () const noexcept { return m_fa.const_arrays(); }

    //! Return mutable multi-array views of every mask.
    MultiArray4<int      > arrays ()       noexcept { return m_fa.arrays(); }

    //! Return const multi-array views (alias of arrays()).
    MultiArray4<int const> const_arrays () const noexcept { return m_fa.const_arrays(); }

    //! Number of components stored in each mask.
    int nComp () const noexcept { return m_fa.nComp(); }

    //! BoxArray used to build this MultiMask.
    const BoxArray& boxArray () const noexcept { return m_fa.boxArray(); }

    //! DistributionMapping describing mask ownership.
    const DistributionMapping& DistributionMap () const noexcept { return m_fa.DistributionMap(); }

    /**
     * \brief Copy multi-mask data from \p src into \p dst.
     *
     * \param dst Destination MultiMask.
     * \param src Source MultiMask.
     */
    static void Copy (MultiMask& dst, const MultiMask& src);

    friend class MultiMaskIter;

private:

    FabArray<Mask> m_fa;
};

//! Iterator adapter exposing MFIter semantics for MultiMask traversal.
class MultiMaskIter
    : public MFIter
{
public:
    explicit MultiMaskIter (const MultiMask& mm)
        : MFIter(mm.m_fa) { }
};

}

#endif
