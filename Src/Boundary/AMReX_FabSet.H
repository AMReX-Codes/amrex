
#ifndef AMREX_FABSET_H_
#define AMREX_FABSET_H_
#include <AMReX_Config.H>

#include <AMReX_FabDataType.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_VisMF.H>

/**
 * \file AMReX_FabSet.H
 * \brief Convenience wrapper for groups of FArrayBoxes located on grid boundaries.
 */

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <limits>

namespace amrex {

/**
        \brief A FabSet is a group of FArrayBox's.  The grouping is designed
        specifically to represent regions along the boundary of Box's,
        and are used to implement boundary conditions to discretized
        partial differential equations.

        A FabSet is an array of pointers to FABs.  The standard FAB operators,
        however, have been modified to be more useful for maintaining
        boundary conditions for partial differential equations discretized
        on boxes.
        Under normal circumstances, a FAB will be created for each face of a
        box.  For a group of boxes, a FabSet will be the group of FABs at a
        particular orientation (ie. the lo-i side of each grid in a list).

        Since a FabSet FAB will likely be used to bound a grid box,
        FArrayBox::resize() operations are disallowed.  Also, to preserve
        flexibility in applicable boundary scenarios, intersecting
        FABs in the FabSet are not guaranteed to contain identical data--thus
        copy operations from a FabSet to any FAB-like structure may be
        order-dependent.

        FabSets are used primarily as a data storage mechanism, and are
        manipulated by more sophisticated control classes.
*/
template <typename MF>
class FabSetT
{
    friend class FabSetIter;
    friend class FluxRegister;
public:
    using value_type = typename FabDataType<MF>::value_type;
    using FAB        = typename FabDataType<MF>::fab_type;

    //
    //! The default constructor -- you must later call define().
    FabSetT () noexcept = default;
    //
    //! Construct a FabSetT<MF> covering \p grids with mapping \p dmap and \p ncomp components.
    FabSetT (const BoxArray& grids, const DistributionMapping& dmap, int ncomp);

    ~FabSetT () = default;

    FabSetT (FabSetT<MF>&& rhs) noexcept = default;

    FabSetT (const FabSetT<MF>& rhs) = delete;
    FabSetT<MF>& operator= (const FabSetT<MF>& rhs) = delete;
    FabSetT<MF>& operator= (FabSetT<MF>&& rhs) = delete;

    //
    //! Define a FabSetT<MF> constructed via default constructor on \p grids/\p dmap with \p ncomp components.
    void define (const BoxArray& grids, const DistributionMapping& dmap, int ncomp);

    //! Access the FAB referenced by iterator \p mfi (const).
    FAB const& operator[] (const MFIter& mfi) const noexcept { return m_mf[mfi]; }

    //! Access the FAB referenced by iterator \p mfi (mutable).
    FAB      & operator[] (const MFIter& mfi)       noexcept { return m_mf[mfi]; }

    //! Access the FAB with index \p i (const).
    FAB const& operator[] (int i)             const noexcept { return m_mf[i]; }

    //! Access the FAB with index \p i (mutable).
    FAB      & operator[] (int i)                   noexcept { return m_mf[i]; }

    //! Return an Array4 view for iterator \p mfi (const).
    Array4<value_type const> array (const MFIter& mfi) const noexcept { return m_mf.const_array(mfi); }

    //! Return an Array4 view for iterator \p mfi (mutable).
    Array4<value_type      > array (const MFIter& mfi)       noexcept { return m_mf.array(mfi); }

    //! Return an Array4 view for FAB index \p i (const).
    Array4<value_type const> array (int i) const noexcept { return m_mf.const_array(i);   }

    //! Return an Array4 view for FAB index \p i (mutable).
    Array4<value_type      > array (int i)       noexcept { return m_mf.array(i);   }

    //! Return a const Array4 view for iterator \p mfi.
    Array4<value_type const> const_array (const MFIter& mfi) const noexcept { return m_mf.const_array(mfi); }

    //! Return a const Array4 view for FAB index \p i.
    Array4<value_type const> const_array (int i)             const noexcept { return m_mf.const_array(i); }

    //! Return multi-fab const Array4 views to every FAB.
    MultiArray4<value_type const> arrays () const noexcept { return m_mf.const_arrays(); }

    //! Return multi-fab mutable Array4 views to every FAB.
    MultiArray4<value_type      > arrays ()       noexcept { return m_mf.arrays(); }

    //! Return multi-fab const Array4 views (alias of arrays()).
    MultiArray4<value_type const> const_arrays () const noexcept { return m_mf.const_arrays(); }

    //! Get the bounding Box of FAB index \p K.
    [[nodiscard]] Box fabbox (int K) const noexcept { return m_mf.fabbox(K); }

    //! Return number of FABs stored in this set.
    [[nodiscard]] int size () const noexcept { return m_mf.size(); }

    //! Return the BoxArray defining FAB locations.
    [[nodiscard]] const BoxArray& boxArray () const noexcept { return m_mf.boxArray(); }

    //! Return the DistributionMapping describing FAB ownership.
    [[nodiscard]] const DistributionMapping& DistributionMap () const noexcept
        { return m_mf.DistributionMap(); }

    //! Mutable access to the underlying MultiFab.
    [[nodiscard]] MF      & multiFab ()       noexcept { return m_mf; }

    //! Const access to the underlying MultiFab.
    [[nodiscard]] MF const& multiFab () const noexcept { return m_mf; }

    //! Number of components carried by each FAB.
    [[nodiscard]] int nComp () const noexcept { return m_mf.nComp(); }

    //! Release storage and reset the internal MultiFab.
    void clear () { m_mf.clear(); }

    /**
     * \brief Copy components from another FabSet with identical layout.
     *
     * \param src   Source FabSet.
     * \param scomp Starting component in \p src.
     * \param dcomp Destination component index.
     * \param ncomp Number of components to copy.
     * \return Reference to this set.
     */
    FabSetT<MF>& copyFrom (const FabSetT<MF>& src, int scomp, int dcomp, int ncomp);

    /**
     * \brief Copy boundary data from a MultiFab intersecting the same boxes.
     *
     * \param src    Source MultiFab.
     * \param ngrow  Number of ghost cells pulled from \p src.
     * \param scomp  Source component index.
     * \param dcomp  Destination component index.
     * \param ncomp  Number of components.
     * \param period Periodicity description.
     */
    FabSetT<MF>& copyFrom (const MF& src, int ngrow, int scomp, int dcomp, int ncomp,
                           const Periodicity& period = Periodicity::NonPeriodic());

    /**
     * \brief Accumulate data from another FabSet with identical layout.
     *
     * \param src   Source FabSet.
     * \param scomp Source component index.
     * \param dcomp Destination component index.
     * \param ncomp Number of components to accumulate.
     */
    FabSetT<MF>& plusFrom (const FabSetT<MF>& src, int scomp, int dcomp, int ncomp);

    /**
     * \brief Accumulate data from a MultiFab into this boundary set.
     *
     * \param src    Source MultiFab.
     * \param ngrow  Number of ghost cells read from \p src.
     * \param scomp  Source component index.
     * \param dcomp  Destination component index.
     * \param ncomp  Number of components.
     * \param period Periodicity description.
     * \param deterministic Whether GPU operations must be deterministic.
     */
    FabSetT<MF>& plusFrom (const MF& src, int ngrow, int scomp, int dcomp, int ncomp,
                           const Periodicity& period = Periodicity::NonPeriodic(),
                           bool deterministic = false);

    /**
     * \brief Copy data from this boundary set back into a MultiFab.
     *
     * \param dest   Destination MultiFab.
     * \param ngrow  Number of grow cells on \p dest to fill.
     * \param scomp  Source component index inside this set.
     * \param dcomp  Destination component index.
     * \param ncomp  Number of components.
     * \param period Periodicity description.
     */
    void copyTo (MF& dest, int ngrow, int scomp, int dcomp, int ncomp,
                 const Periodicity& period = Periodicity::NonPeriodic()) const;

    /**
     * \brief Add this boundary data back into a MultiFab.
     *
     * \param dest   Destination MultiFab.
     * \param ngrow  Number of grow cells to cover on \p dest.
     * \param scomp  Source component index inside this set.
     * \param dcomp  Destination component index.
     * \param ncomp  Number of components.
     * \param period Periodicity description.
     */
    void plusTo (MF& dest, int ngrow, int scomp, int dcomp, int ncomp,
                 const Periodicity& period = Periodicity::NonPeriodic()) const;

    //! Fill every component of every FAB with scalar value \p val.
    void setVal (value_type val);

    //! Fill \p num_comp components starting at \p comp with scalar value \p val.
    void setVal (value_type val, int comp, int num_comp);

    /**
     * \brief Linear combination: this := a*this + b*src.
     *
     * \param a     Scaling applied to the current FabSet.
     * \param b     Scaling applied to \p src.
     * \param src   Source FabSet (must match layout).
     * \param scomp Source component offset.
     * \param dcomp Destination component offset.
     * \param ncomp Number of components.
     */
    FabSetT<MF>& linComb (value_type a, value_type b, const FabSetT<MF>& src,
                          int scomp, int dcomp, int ncomp);

    /**
     * \brief Boundary linear combination from two MultiFabs: this := a*mfa + b*mfb.
     *
     * \param a      Scale applied to \p mfa.
     * \param mfa    First MultiFab input.
     * \param a_comp Component offset in \p mfa.
     * \param b      Scale applied to \p mfb.
     * \param mfb    Second MultiFab input.
     * \param b_comp Component offset in \p mfb.
     * \param dcomp  Destination component offset.
     * \param ncomp  Number of components.
     * \param ngrow  Number of grow cells required in inputs.
     */
    FabSetT<MF>& linComb (value_type a, const MF& mfa, int a_comp,
                          value_type b, const MF& mfb, int b_comp,
                          int dcomp, int ncomp, int ngrow);

    /**
     * \brief Write (used for writing to checkpoint).
     *
     * \param name Base filename to write.
     */
    void write (const std::string& name) const;

    /**
     * \brief Read (used for reading from checkpoint).
     *
     * \param name Base filename previously written by write().
     */
    void read (const std::string& name);

    /**
     * \brief Local copy function.
     *
     * \param dst Destination FabSet.
     * \param src Source FabSet.
     */
    static void Copy (FabSetT<MF>& dst, const FabSetT<MF>& src);

private:
    MF m_mf;
};

//! Convenience iterator that reuses MFIter semantics for FabSet traversal.
class FabSetIter
    : public MFIter
{
public:
    template <typename MF>
    explicit FabSetIter (const FabSetT<MF>& fs)
        : MFIter(fs.m_mf) { }
};

template <typename MF>
FabSetT<MF>::FabSetT (const BoxArray& grids, const DistributionMapping& dmap, int ncomp)
    : m_mf(grids,dmap,ncomp,0)
{}

template <typename MF>
void
FabSetT<MF>::define (const BoxArray& grids, const DistributionMapping& dm, int ncomp)
{
    m_mf.define(grids, dm, ncomp, 0);
}

template <typename MF>
FabSetT<MF>&
FabSetT<MF>::copyFrom (const FabSetT<MF>& src, int scomp, int dcomp, int ncomp)
{
    if (boxArray() == src.boxArray() && DistributionMap() == src.DistributionMap()) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (FabSetIter fsi(*this); fsi.isValid(); ++fsi) {
            const Box& bx = fsi.validbox();
            auto const srcfab =   src.array(fsi);
            auto       dstfab = this->array(fsi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                dstfab(i,j,k,n+dcomp) = srcfab(i,j,k,n+scomp);
            });
        }
    } else {
        m_mf.ParallelCopy(src.m_mf,scomp,dcomp,ncomp);
    }
    return *this;
}

template <typename MF>
FabSetT<MF>&
FabSetT<MF>::copyFrom (const MF& src, int ngrow, int scomp, int dcomp, int ncomp,
                       const Periodicity& period)
{
    BL_ASSERT(boxArray() != src.boxArray());
    m_mf.ParallelCopy(src,scomp,dcomp,ncomp,ngrow,0,period);
    return *this;
}

template <typename MF>
FabSetT<MF>&
FabSetT<MF>::plusFrom (const FabSetT<MF>& src, int scomp, int dcomp, int ncomp)
{
    if (boxArray() == src.boxArray() && DistributionMap() == src.DistributionMap()) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (FabSetIter fsi(*this); fsi.isValid(); ++fsi) {
            const Box& bx = fsi.validbox();
            auto const srcfab =   src.array(fsi);
            auto       dstfab = this->array(fsi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                dstfab(i,j,k,n+dcomp) += srcfab(i,j,k,n+scomp);
            });
        }
    } else {
        amrex::Abort("FabSetT<MF>::plusFrom: parallel plusFrom not supported");
    }
    return *this;
}

template <typename MF>
FabSetT<MF>&
FabSetT<MF>::plusFrom (const MF& src, int ngrow, int scomp, int dcomp, int ncomp,
                       const Periodicity& period, bool deterministic)
{
    BL_ASSERT(boxArray() != src.boxArray());
    m_mf.ParallelCopy(src,scomp,dcomp,ncomp,IntVect(ngrow),IntVect(0),period,
                      FabArrayBase::ADD,nullptr,deterministic);
    return *this;
}

template <typename MF>
void
FabSetT<MF>::copyTo (MF& dest, int ngrow, int scomp, int dcomp, int ncomp,
                     const Periodicity& period) const
{
    BL_ASSERT(boxArray() != dest.boxArray());
    dest.ParallelCopy(m_mf,scomp,dcomp,ncomp,0,ngrow,period);
}

template <typename MF>
void
FabSetT<MF>::plusTo (MF& dest, int ngrow, int scomp, int dcomp, int ncomp,
                     const Periodicity& period) const
{
    BL_ASSERT(boxArray() != dest.boxArray());
    dest.ParallelCopy(m_mf,scomp,dcomp,ncomp,0,ngrow,period,FabArrayBase::ADD);
}

template <typename MF>
void
FabSetT<MF>::setVal (value_type val)
{
    const int ncomp = nComp();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi) {
        const Box& bx = fsi.validbox();
        auto fab = this->array(fsi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            fab(i,j,k,n) = val;
        });
    }
}

template <typename MF>
void
FabSetT<MF>::setVal (value_type val, int comp, int num_comp)
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi) {
        const Box& bx = fsi.validbox();
        auto fab = this->array(fsi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, num_comp, i, j, k, n,
        {
            fab(i,j,k,n+comp) = val;
        });
    }
}

// Linear combination this := a*this + b*src
// Note: corresponding fabsets must be commensurate.
template <typename MF>
FabSetT<MF>&
FabSetT<MF>::linComb (value_type a, value_type b, const FabSetT<MF>& src,
                      int scomp, int dcomp, int ncomp)
{
    BL_ASSERT(size() == src.size());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi)
    {
        const Box& bx = fsi.validbox();
        auto const srcfab =   src.array(fsi);
        auto       dstfab = this->array(fsi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            dstfab(i,j,k,n+dcomp) = a*dstfab(i,j,k,n+dcomp) + b*srcfab(i,j,k,n+scomp);
        });
    }
    return *this;
}

// Linear combination: this := a*mfa + b*mfb
// CastroRadiation is the only code that uses this function.
template <typename MF>
FabSetT<MF>&
FabSetT<MF>::linComb (value_type a, const MF& mfa, int a_comp,
                      value_type b, const MF& mfb, int b_comp,
                      int dcomp, int ncomp, int ngrow)
{
    BL_PROFILE("FabSetT<MF>::linComb()");
    BL_ASSERT(mfa.nGrowVect().allGE(ngrow) && mfb.nGrowVect().allGE(ngrow));
    BL_ASSERT(mfa.boxArray() == mfb.boxArray());
    BL_ASSERT(boxArray() != mfa.boxArray());

    MF bdrya(boxArray(),DistributionMap(),ncomp,0,MFInfo());
    MF bdryb(boxArray(),DistributionMap(),ncomp,0,MFInfo());

    const auto huge = static_cast<value_type>(sizeof(value_type) == 8 ? 1.e200 : 1.e30);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(bdrya); mfi.isValid(); ++mfi) // tiling is not safe for this BoxArray
    {
        const Box& bx = mfi.validbox();
        auto afab = bdrya.array(mfi);
        auto bfab = bdryb.array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            afab(i,j,k,n) = huge;
            bfab(i,j,k,n) = huge;
        });
    }

    bdrya.ParallelCopy(mfa,a_comp,0,ncomp,ngrow,0);
    bdryb.ParallelCopy(mfb,b_comp,0,ncomp,ngrow,0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi)
    {
        const Box& bx = fsi.validbox();
        auto const afab = bdrya.array(fsi);
        auto const bfab = bdryb.array(fsi);
        auto       dfab = this->array(fsi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            dfab(i,j,k,n+dcomp) = a*afab(i,j,k,n) + b*bfab(i,j,k,n);
        });
    }

    return *this;
}

template <typename MF>
void
FabSetT<MF>::write (const std::string& name) const
{
    if (AsyncOut::UseAsyncOut()) {
        VisMF::AsyncWrite(m_mf,name);
    } else {
        VisMF::Write(m_mf,name);
    }
}

template <typename MF>
void
FabSetT<MF>::read (const std::string& name)
{
    if (m_mf.empty()) {
        amrex::Abort("FabSetT<MF>::read: not predefined");
    }
    VisMF::Read(m_mf,name);
}

template <typename MF>
void
FabSetT<MF>::Copy (FabSetT<MF>& dst, const FabSetT<MF>& src)
{
    BL_ASSERT(amrex::match(dst.boxArray(), src.boxArray()));
    BL_ASSERT(dst.DistributionMap() == src.DistributionMap());
    int ncomp = dst.nComp();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (FabSetIter fsi(dst); fsi.isValid(); ++fsi) {
        const Box& bx = fsi.validbox();
        auto const srcfab = src.array(fsi);
        auto       dstfab = dst.array(fsi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            dstfab(i,j,k,n) = srcfab(i,j,k,n);
        });
    }
}

using FabSet = FabSetT<MultiFab>;
using fFabSet = FabSetT<fMultiFab>;

}

#endif /*_FABSET_H_*/
