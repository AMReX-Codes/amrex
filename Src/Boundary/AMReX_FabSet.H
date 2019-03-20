
#ifndef AMREX_FABSET_H_
#define AMREX_FABSET_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

namespace amrex {

class FabSetIter;

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

class FabSet
{
    friend class FabSetIter;
    friend class FluxRegister;
public:
    //
    //! The default constructor -- you must later call define().
    FabSet () noexcept;
    //
    //! Construct a FabSet of specified number of components on the grids.
    FabSet (const BoxArray& grids, const DistributionMapping& dmap, int ncomp);

    ~FabSet () = default;

    FabSet (FabSet&& rhs) noexcept = default;

    FabSet (const FabSet& rhs) = delete;
    FabSet& operator= (const FabSet& rhs) = delete;
    FabSet& operator= (FabSet&& rhs) = delete;

    //
    //! Define a FabSet constructed via default constructor.
    void define (const BoxArray& grids, const DistributionMapping& dmap, int ncomp);

    FArrayBox const& operator[] (const MFIter& mfi) const noexcept { return m_mf[mfi]; }
    FArrayBox      & operator[] (const MFIter& mfi)       noexcept { return m_mf[mfi]; }
    FArrayBox const& operator[] (int i)             const noexcept { return m_mf[i]; }
    FArrayBox      & operator[] (int i)                   noexcept { return m_mf[i]; }

    FArrayBox const* fabPtr (const MFIter& mfi) const noexcept { return m_mf.fabPtr(mfi); }
    FArrayBox      * fabPtr (const MFIter& mfi)       noexcept { return m_mf.fabPtr(mfi); }
    FArrayBox const* fabPtr (int i)             const noexcept { return m_mf.fabPtr(i); }
    FArrayBox      * fabPtr (int i)                   noexcept { return m_mf.fabPtr(i); }

    FArrayBox const* fabHostPtr (const MFIter& mfi) const noexcept { return m_mf.fabHostPtr(mfi); }
    FArrayBox      * fabHostPtr (const MFIter& mfi)       noexcept { return m_mf.fabHostPtr(mfi); }
    FArrayBox const* fabHostPtr (int i)             const noexcept { return m_mf.fabHostPtr(i); }
    FArrayBox      * fabHostPtr (int i)                   noexcept { return m_mf.fabHostPtr(i); }

    FArrayBox const* fabDevicePtr (const MFIter& mfi) const noexcept { return m_mf.fabDevicePtr(mfi); }
    FArrayBox      * fabDevicePtr (const MFIter& mfi)       noexcept { return m_mf.fabDevicePtr(mfi); }
    FArrayBox const* fabDevicePtr (int i)             const noexcept { return m_mf.fabDevicePtr(i); }
    FArrayBox      * fabDevicePtr (int i)                   noexcept { return m_mf.fabDevicePtr(i); }

    Array4<Real const> array (const MFIter& mfi) const noexcept { return m_mf.array(mfi); }
    Array4<Real      > array (const MFIter& mfi)       noexcept { return m_mf.array(mfi); }
    Array4<Real const> array (int i)             const noexcept { return m_mf.array(i);   }
    Array4<Real      > array (int i)                   noexcept { return m_mf.array(i);   }

    Box fabbox (int K) const noexcept { return m_mf.fabbox(K); }

    int size () const noexcept { return m_mf.size(); }

    const BoxArray& boxArray () const noexcept { return m_mf.boxArray(); }

    const DistributionMapping& DistributionMap () const noexcept
	{ return m_mf.DistributionMap(); }

    int nComp () const noexcept { return m_mf.nComp(); }

    void clear () { m_mf.clear(); }

    FabSet& copyFrom (const FabSet& src, int scomp, int dcomp, int ncomp);

    FabSet& copyFrom (const MultiFab& src, int ngrow, int scomp, int dcomp, int ncomp,
                      const Periodicity& period = Periodicity::NonPeriodic());

    FabSet& plusFrom (const FabSet& src, int scomp, int dcomp, int ncomp);

    FabSet& plusFrom (const MultiFab& src, int ngrow, int scomp, int dcomp, int ncomp,
		      const Periodicity& period = Periodicity::NonPeriodic());

    void copyTo (MultiFab& dest, int ngrow, int scomp, int dcomp, int ncomp,
		 const Periodicity& period = Periodicity::NonPeriodic()) const;

    void plusTo (MultiFab& dest, int ngrow, int scomp, int dcomp, int ncomp,
		 const Periodicity& period = Periodicity::NonPeriodic()) const;

    void setVal (Real val);

    void setVal (Real val, int comp, int num_comp);

    //!< Linear combination: this := a*this + b*src (FabSets must be commensurate).
    FabSet& linComb (Real a, Real b, const FabSet& src, int scomp, int dcomp, int ncomp);

    //!< Linear combination: this := a*mfa + b*mfb
    FabSet& linComb (Real a, const MultiFab& mfa, int a_comp,
		     Real b, const MultiFab& mfb, int b_comp,
                     int dcomp, int ncomp, int ngrow);

    //
    //! Write (used for writing to checkpoint)
    void write (const std::string& name) const;
    //
    //! Read (used for reading from checkpoint)
    void read (const std::string& name);

    //!< Local copy function
    static void Copy (FabSet& dst, const FabSet& src);

private:
    MultiFab m_mf;
};

class FabSetIter
    : public MFIter
{
public:
    explicit FabSetIter (const FabSet& fs)
	: MFIter(fs.m_mf) { }
};

}

#endif /*_FABSET_H_*/
