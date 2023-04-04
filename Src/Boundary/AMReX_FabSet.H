
#ifndef AMREX_FABSET_H_
#define AMREX_FABSET_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_VisMF.H>

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
    using value_type = typename MF::value_type;
    using FAB = typename MF::fab_type;

    //
    //! The default constructor -- you must later call define().
    FabSetT () noexcept = default;
    //
    //! Construct a FabSetT<MF> of specified number of components on the grids.
    FabSetT (const BoxArray& grids, const DistributionMapping& dmap, int ncomp);

    ~FabSetT () = default;

    FabSetT (FabSetT<MF>&& rhs) noexcept = default;

    FabSetT (const FabSetT<MF>& rhs) = delete;
    FabSetT<MF>& operator= (const FabSetT<MF>& rhs) = delete;
    FabSetT<MF>& operator= (FabSetT<MF>&& rhs) = delete;

    //
    //! Define a FabSetT<MF> constructed via default constructor.
    void define (const BoxArray& grids, const DistributionMapping& dmap, int ncomp);

    FAB const& operator[] (const MFIter& mfi) const noexcept { return m_mf[mfi]; }
    FAB      & operator[] (const MFIter& mfi)       noexcept { return m_mf[mfi]; }
    FAB const& operator[] (int i)             const noexcept { return m_mf[i]; }
    FAB      & operator[] (int i)                   noexcept { return m_mf[i]; }

    Array4<value_type const> array (const MFIter& mfi) const noexcept { return m_mf.const_array(mfi); }
    Array4<value_type      > array (const MFIter& mfi)       noexcept { return m_mf.array(mfi); }
    Array4<value_type const> array (int i) const noexcept { return m_mf.const_array(i);   }
    Array4<value_type      > array (int i)       noexcept { return m_mf.array(i);   }
    Array4<value_type const> const_array (const MFIter& mfi) const noexcept { return m_mf.const_array(mfi); }
    Array4<value_type const> const_array (int i)             const noexcept { return m_mf.const_array(i); }

    MultiArray4<value_type const> arrays () const noexcept { return m_mf.const_arrays(); }
    MultiArray4<value_type      > arrays ()       noexcept { return m_mf.arrays(); }
    MultiArray4<value_type const> const_arrays () const noexcept { return m_mf.const_arrays(); }

    Box fabbox (int K) const noexcept { return m_mf.fabbox(K); }

    int size () const noexcept { return m_mf.size(); }

    const BoxArray& boxArray () const noexcept { return m_mf.boxArray(); }

    const DistributionMapping& DistributionMap () const noexcept
        { return m_mf.DistributionMap(); }

    MF      & multiFab ()       noexcept { return m_mf; }
    MF const& multiFab () const noexcept { return m_mf; }

    int nComp () const noexcept { return m_mf.nComp(); }

    void clear () { m_mf.clear(); }

    FabSetT<MF>& copyFrom (const FabSetT<MF>& src, int scomp, int dcomp, int ncomp);

    FabSetT<MF>& copyFrom (const MF& src, int ngrow, int scomp, int dcomp, int ncomp,
                           const Periodicity& period = Periodicity::NonPeriodic());

    FabSetT<MF>& plusFrom (const FabSetT<MF>& src, int scomp, int dcomp, int ncomp);

    FabSetT<MF>& plusFrom (const MF& src, int ngrow, int scomp, int dcomp, int ncomp,
                           const Periodicity& period = Periodicity::NonPeriodic());

    void copyTo (MF& dest, int ngrow, int scomp, int dcomp, int ncomp,
                 const Periodicity& period = Periodicity::NonPeriodic()) const;

    void plusTo (MF& dest, int ngrow, int scomp, int dcomp, int ncomp,
                 const Periodicity& period = Periodicity::NonPeriodic()) const;

    void setVal (value_type val);

    void setVal (value_type val, int comp, int num_comp);

    //!< Linear combination: this := a*this + b*src (FabSetT<MF>s must be commensurate).
    FabSetT<MF>& linComb (value_type a, value_type b, const FabSetT<MF>& src,
                          int scomp, int dcomp, int ncomp);

    //!< Linear combination: this := a*mfa + b*mfb
    FabSetT<MF>& linComb (value_type a, const MF& mfa, int a_comp,
                          value_type b, const MF& mfb, int b_comp,
                          int dcomp, int ncomp, int ngrow);

    //
    //! Write (used for writing to checkpoint)
    void write (const std::string& name) const;
    //
    //! Read (used for reading from checkpoint)
    void read (const std::string& name);

    //!< Local copy function
    static void Copy (FabSetT<MF>& dst, const FabSetT<MF>& src);

private:
    MF m_mf;
};

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
                       const Periodicity& period)
{
    BL_ASSERT(boxArray() != src.boxArray());
    m_mf.ParallelCopy(src,scomp,dcomp,ncomp,ngrow,0,period,FabArrayBase::ADD);
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
    BL_ASSERT(ngrow <= mfa.nGrow());
    BL_ASSERT(ngrow <= mfb.nGrow());
    BL_ASSERT(mfa.boxArray() == mfb.boxArray());
    BL_ASSERT(boxArray() != mfa.boxArray());

    MF bdrya(boxArray(),DistributionMap(),ncomp,0,MFInfo());
    MF bdryb(boxArray(),DistributionMap(),ncomp,0,MFInfo());

    const value_type huge = static_cast<value_type>(sizeof(value_type) == 8 ? 1.e200 : 1.e30);

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
