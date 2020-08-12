
#ifndef BL_FABARRAY_H
#define BL_FABARRAY_H

#include <iostream>
#include <cstring>
#include <limits>
#include <map>
#include <utility>
#include <vector>
#include <algorithm>
#include <set>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_BLassert.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_Box.H>
#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_BoxDomain.H>
#include <AMReX_FabFactory.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_ccse-mpi.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Periodicity.H>
#include <AMReX_Print.H>
#include <AMReX_FabArrayBase.H>
#include <AMReX_MFIter.H>
#include <AMReX_MakeType.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_LayoutData.H>
#include <AMReX_BaseFab.H>

#include <AMReX_Gpu.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif

namespace amrex {

template <typename T, typename std::enable_if<!IsBaseFab<T>::value,int>::type = 0>
Long nBytesOwned (T const&) noexcept { return 0; }

template <typename T>
Long nBytesOwned (BaseFab<T> const& fab) noexcept { return fab.nBytesOwned(); }

/*
  A Collection of Fortran Array-like Objects


  The FabArray<FAB> class implements a collection (stored as an array) of
  Fortran array-like objects.  The parameterized type FAB is intended to be
  any class derived from BaseFab<T>.  For example, FAB may be a BaseFab of
  integers, so we could write:

    FabArray<BaseFab<int> > int_fabs;

  Then int_fabs is a FabArray that can hold a collection of BaseFab<int>
  objects.

  FabArray is not just a general container class for Fortran arrays.  It is
  intended to hold "grid" data for use in finite difference calculations in
  which the data is defined on a union of (usually disjoint) rectangular
  regions embedded in a uniform index space.  This region, called the valid
  region, is represented by a BoxArray.  For the purposes of this discussion,
  the Kth Box in the BoxArray represents the interior region of the Kth grid.

  Since the intent is to be used with finite difference calculations a
  FabArray also includes the notion of a boundary region for each grid.  The
  boundary region is specified by the ngrow parameter which tells the FabArray
  to allocate each FAB to be ngrow cells larger in all directions than the
  underlying Box.  The larger region covered by the union of all the FABs is
  called the region of definition.  The underlying notion is that the valid
  region contains the grid interior data and the region of definition includes
  the interior region plus the boundary areas.

  Operations are available to copy data from the valid regions into these
  boundary areas where the two overlap.  The number of components, that is,
  the number of values that can be stored in each cell of a FAB, is either
  given as an argument to the constructor or is inherent in the definition of
  the underlying FAB.  Each FAB in the FabArray will have the same number of
  components.

  In summary, a FabArray is an array of FABs.  The Kth element contains a FAB
  that holds the data for the Kth grid, a Box that defines the valid region
  of the Kth grid.

  A typical use for a FabArray would be to hold the solution vector or
  right-hand-side when solving a linear system of equations on a union of
  rectangular grids.  The copy operations would be used to copy data from the
  valid regions of neighboring grids into the boundary regions after each
  relaxation step of the iterative method.  If a multigrid method is used, a
  FabArray could be used to hold the data at each level in the multigrid
  hierarchy.

  This class is a concrete class not a polymorphic one.

  This class does NOT provide a copy constructor or assignment operator.
*/

//
// alloc: allocate memory or not
//
struct MFInfo {
    bool    alloc = true;
    Arena*  arena = nullptr;
    Vector<std::string> tags;

    MFInfo& SetAlloc (bool a) noexcept { alloc = a; return *this; }

    MFInfo& SetArena (Arena* ar) noexcept { arena = ar; return *this; }

    MFInfo& SetTag (const char* t) noexcept {
        tags.emplace_back(t);
        return *this;
    }

    MFInfo& SetTag (std::string t) noexcept {
        tags.emplace_back(std::move(t));
        return *this;
    }

    template <typename T, typename... Ts>
    MFInfo& SetTag (T&& t, Ts&&... ts) noexcept {
        tags.emplace_back(std::forward<T>(t));
        return SetTag(std::forward<Ts>(ts)...);
    }
};

template <class FAB>
class FabArray
    :
    public FabArrayBase
{
public:

    struct FABType {
        typedef FAB value_type;
    };

    /*
    * if FAB is a BaseFab or its child, value_type = FAB::value_type
    * else                              value_type = FAB;
    */
    using value_type = typename std::conditional<IsBaseFab<FAB>::value, FAB, FABType>::type::value_type;

    //
    //! Constructs an empty FabArray<FAB>.
    FabArray ();

    /**
    * \brief Construct a FabArray<FAB> with a valid region defined by bxs
    * and a region of definition defined by the grow factor ngrow
    * and the number of components nvar.
    */
    FabArray (const BoxArray&            bxs,
              const DistributionMapping& dm,
              int                        nvar,
              int                        ngrow,
#ifdef AMREX_STRICT_MODE
	      const MFInfo&              info,
              const FabFactory<FAB>&     factory);
#else
	      const MFInfo&              info = MFInfo(),
              const FabFactory<FAB>&     factory = DefaultFabFactory<FAB>());
#endif

    FabArray (const BoxArray&            bxs,
              const DistributionMapping& dm,
              int                        nvar,
              const IntVect&             ngrow,
#ifdef AMREX_STRICT_MODE
	      const MFInfo&              info,
              const FabFactory<FAB>&     factory);
#else
	      const MFInfo&              info = MFInfo(),
              const FabFactory<FAB>&     factory = DefaultFabFactory<FAB>());
#endif

    FabArray (const FabArray<FAB>& rhs, MakeType maketype, int scomp, int ncomp);

    //! The destructor -- deletes all FABs in the array.
    virtual ~FabArray ();

    FabArray (FabArray<FAB>&& rhs) noexcept;
    FabArray<FAB>& operator= (FabArray<FAB>&& rhs) noexcept;

    FabArray (const FabArray<FAB>& rhs) = delete;
    FabArray<FAB>& operator= (const FabArray<FAB>& rhs) = delete;

    /**
    * \brief Define this FabArray identically to that performed by
    * the constructor having an analogous function signature.
    * This is only valid if this FabArray was defined using
    * the default constructor.
    */
    virtual void define (const BoxArray& bxs,
			 const DistributionMapping& dm,
			 int                        nvar,
			 int                        ngrow,
#ifdef AMREX_STRICT_MODE
			 const MFInfo&              info,
                         const FabFactory<FAB>&     factory);
#else
			 const MFInfo&              info = MFInfo(),
                         const FabFactory<FAB>&     factory = DefaultFabFactory<FAB>());
#endif

    virtual void define (const BoxArray& bxs,
			 const DistributionMapping& dm,
			 int                        nvar,
			 const IntVect&             ngrow,
#ifdef AMREX_STRICT_MODE
			 const MFInfo&              info,
                         const FabFactory<FAB>&     factory);
#else
			 const MFInfo&              info = MFInfo(),
                         const FabFactory<FAB>&     factory = DefaultFabFactory<FAB>());
#endif

    const FabFactory<FAB>& Factory () const noexcept { return *m_factory; }

    // Provides access to the Arena this FabArray was build with.
    Arena* arena () const noexcept { return m_dallocator.arena(); }

    bool hasEBFabFactory () const noexcept {
#ifdef AMREX_USE_EB
        const auto f = dynamic_cast<EBFArrayBoxFactory const*>(m_factory.get());
        return (f != nullptr);
#else
        return false;
#endif
    }

    bool isAllRegular () const noexcept {
#ifdef AMREX_USE_EB
        const auto f = dynamic_cast<EBFArrayBoxFactory const*>(m_factory.get());
        if (f) {
            return f->isAllRegular();
        } else {
            return true;
        }
#else
        return true;
#endif
    }

    /**
    * \brief Return true if the FabArray is well-defined.  That is,
    * the FabArray has a BoxArray and DistributionMapping, the
    * FABs are allocated for each Box in the BoxArray and the
    * sizes of the FABs and the number of components are consistent
    * with the definition of the FabArray.
    */
    bool ok () const;

    //! Return a constant reference to the FAB associated with mfi.
    const FAB& operator[] (const MFIter& mfi) const noexcept { return *(this->fabPtr(mfi)); }

    //! Return a constant reference to the FAB associated with mfi.
    const FAB& get (const MFIter& mfi) const noexcept { return *(this->fabPtr(mfi)); }

    //! Returns a reference to the FAB associated mfi.
    FAB& operator[] (const MFIter& mfi) noexcept { return *(this->fabPtr(mfi)); }

    //! Returns a reference to the FAB associated mfi.
    FAB& get (const MFIter& mfi) noexcept { return *(this->fabPtr(mfi)); }

    //! Return a constant reference to the FAB associated with the Kth element.
    const FAB& operator[] (int K) const noexcept { return *(this->fabPtr(K)); }

    //! Return a constant reference to the FAB associated with the Kth element.
    const FAB& get (int K) const noexcept { return *(this->fabPtr(K)); }

    //! Return a reference to the FAB associated with the Kth element.
    FAB& operator[] (int K) noexcept { return *(this->fabPtr(K)); }

    //! Return a reference to the FAB associated with the Kth element.
    FAB& get (int K) noexcept { return *(this->fabPtr(K)); }

    //! Return a reference to the FAB associated with local index L
    FAB& atLocalIdx (int L) noexcept { return *m_fabs_v[L]; }
    const FAB& atLocalIdx (int L) const noexcept { return *m_fabs_v[L]; }

    //! Return pointer to FAB
    FAB      * fabPtr (const MFIter& mfi) noexcept;
    FAB const* fabPtr (const MFIter& mfi) const noexcept;
    FAB      * fabPtr (int K) noexcept;  // Here K is global index
    FAB const* fabPtr (int K) const noexcept;

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void prefetchToHost (const MFIter& mfi) const noexcept;

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void prefetchToDevice (const MFIter& mfi) const noexcept;

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type const> array (const MFIter& mfi) const noexcept;
    //
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type> array (const MFIter& mfi) noexcept;
    //
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type const> array (int K) const noexcept;
    //
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type> array (int K) noexcept;

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type const> const_array (const MFIter& mfi) const noexcept;
    //
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type const> const_array (int K) const noexcept;

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type const> array (const MFIter& mfi, int start_comp) const noexcept;
    //
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type> array (const MFIter& mfi, int start_comp) noexcept;
    //
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type const> array (int K, int start_comp) const noexcept;
    //
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type> array (int K, int start_comp) noexcept;

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type const> const_array (const MFIter& mfi, int start_comp) const noexcept;
    //
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    Array4<typename FabArray<FAB>::value_type const> const_array (int K, int start_comp) const noexcept;

    //! Explicitly set the Kth FAB in the FabArray to point to elem.
    void setFab (int K, FAB* elem);

    //! Explicitly set the FAB associated with mfi in the FabArray to point to elem.
    void setFab (const MFIter&mfi, FAB* elem, bool assertion=true);

    //! Releases FAB memory in the FabArray.
    void clear ();

    //! Set all components in the entire region of each FAB to val.
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val);

    //! Set all components in the entire region of each FAB to val.
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void operator= (value_type val);

    /**
    * \brief Set the value of num_comp components in the valid region of
    * each FAB in the FabArray, starting at component comp to val.
    * Also set the value of nghost boundary cells.
    */
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val,
                 int        comp,
                 int        num_comp,
                 const IntVect& nghost);

    /**
    * \brief Set the value of num_comp components in the valid region of
    * each FAB in the FabArray, starting at component comp, as well
    * as nghost boundary cells, to val, provided they also intersect
    * with the Box region.
    */
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val,
                 const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val,
                 const Box& region,
                 int        comp,
                 int        num_comp,
                 const IntVect& nghost);
    /**
    * \brief Set all components in the valid region of each FAB in the
    * FabArray to val, including nghost boundary cells.
    */
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val, int nghost);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val, const IntVect& nghost);

    /**
    * \brief Set all components in the valid region of each FAB in the
    * FabArray to val, including nghost boundary cells, that also
    * intersect the Box region.
    */
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val, const Box& region, int nghost);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type val, const Box& region, const IntVect& nghost);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void abs (int comp, int num_comp, int nghost = 0);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void abs (int comp, int num_comp, const IntVect& nghost);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void plus (value_type val, int comp, int num_comp, int nghost = 0);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void plus (value_type val, const Box& region, int comp, int num_comp, int nghost = 0);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void mult (value_type val, int comp, int num_comp, int nghost = 0);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void mult (value_type val, const Box& region, int comp, int num_comp, int nghost = 0);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void invert (value_type numerator, int comp, int num_comp, int nghost = 0);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void invert (value_type numerator, const Box& region, int comp, int num_comp, int nghost = 0);

    //! Set all values in the boundary region to val.
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setBndry (value_type val);

    //! Set ncomp values in the boundary region, starting at start_comp to val.
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setBndry (value_type val, int strt_comp, int ncomp);

   //! Set all values outside the Geometry domain to val.
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setDomainBndry (value_type val, const Geometry& goem);

    //! Set ncomp values outside the Geometry domain to val, starting at start_comp.
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setDomainBndry (value_type val, int strt_comp, int ncomp, const Geometry& goem);

    /**
    * \brief This function copies data from fa to this FabArray.  Each FAB
    * in fa is intersected with all FABs in this FabArray and a copy
    * is performed on the region of intersection.  The intersection
    * is restricted to the valid regions.
    */
    void ParallelAdd (const FabArray<FAB>& fa,
                      const Periodicity&   period = Periodicity::NonPeriodic())
       { ParallelCopy(fa,period,FabArray::ADD); }
    void ParallelCopy (const FabArray<FAB>& fa,
                       const Periodicity&   period = Periodicity::NonPeriodic(),
                       CpOp                 op = FabArrayBase::COPY)
       { ParallelCopy(fa,0,0,nComp(),0,0,period,op); }
    void copy (const FabArray<FAB>& fa,
	       const Periodicity&   period = Periodicity::NonPeriodic(),
               CpOp                 op = FabArrayBase::COPY)
        { ParallelCopy(fa,period,op); }

    /**
    * \brief This function copies data from src to this FabArray.  Each FAB
    * in src is intersected with all FABs in this FabArray and a copy
    * is performed on the region of intersection.  The intersection
    * is restricted to the num_comp components starting at src_comp
    * in the FabArray src, with the destination components in this
    * FabArray starting at dest_comp.
    */
    void ParallelAdd (const FabArray<FAB>& src,
                      int                  src_comp,
                      int                  dest_comp,
                      int                  num_comp,
                      const Periodicity&   period = Periodicity::NonPeriodic())
       { ParallelCopy(src,src_comp,dest_comp,num_comp, period, FabArrayBase::ADD); }
    void ParallelCopy (const FabArray<FAB>& src,
                       int                  src_comp,
                       int                  dest_comp,
                       int                  num_comp,
                       const Periodicity&   period = Periodicity::NonPeriodic(),
                       CpOp                 op = FabArrayBase::COPY)
       { ParallelCopy(src,src_comp,dest_comp,num_comp,0,0,period,op); }
    void copy (const FabArray<FAB>& src,
               int                  src_comp,
               int                  dest_comp,
               int                  num_comp,
	       const Periodicity&   period = Periodicity::NonPeriodic(),
               CpOp                 op = FabArrayBase::COPY)
        { ParallelCopy(src,src_comp,dest_comp,num_comp, period, op); }

    //! Similar to the above function, except that source and destination are grown by src_nghost and dst_nghost, respectively
    void ParallelAdd (const FabArray<FAB>& src,
                      int                  src_comp,
                      int                  dest_comp,
                      int                  num_comp,
                      int                  src_nghost,
                      int                  dst_nghost,
                      const Periodicity&   period = Periodicity::NonPeriodic())
       { ParallelCopy(src,src_comp,dest_comp,num_comp,IntVect(src_nghost),IntVect(dst_nghost),period,
                      FabArrayBase::ADD); }
    void ParallelAdd (const FabArray<FAB>& src,
                      int                  src_comp,
                      int                  dest_comp,
                      int                  num_comp,
                      const IntVect&       src_nghost,
                      const IntVect&       dst_nghost,
                      const Periodicity&   period = Periodicity::NonPeriodic())
       { ParallelCopy(src,src_comp,dest_comp,num_comp,src_nghost,dst_nghost,period,FabArrayBase::ADD); }
    void ParallelCopy (const FabArray<FAB>& src,
                       int                  src_comp,
                       int                  dest_comp,
                       int                  num_comp,
                       int                  src_nghost,
                       int                  dst_nghost,
                       const Periodicity&   period = Periodicity::NonPeriodic(),
                       CpOp                 op = FabArrayBase::COPY)
       { ParallelCopy(src,src_comp,dest_comp,num_comp,IntVect(src_nghost),IntVect(dst_nghost),period,op); }
    void ParallelCopy (const FabArray<FAB>& src,
                       int                  src_comp,
                       int                  dest_comp,
                       int                  num_comp,
                       const IntVect&       src_nghost,
                       const IntVect&       dst_nghost,
                       const Periodicity&   period = Periodicity::NonPeriodic(),
                       CpOp                 op = FabArrayBase::COPY,
                       const FabArrayBase::CPC* a_cpc = nullptr);

    void copy (const FabArray<FAB>& src,
               int                  src_comp,
               int                  dest_comp,
               int                  num_comp,
	       int                  src_nghost,
	       int                  dst_nghost,
	       const Periodicity&   period = Periodicity::NonPeriodic(),
               CpOp                 op = FabArrayBase::COPY)
       { ParallelCopy(src,src_comp,dest_comp,num_comp,IntVect(src_nghost),IntVect(dst_nghost),period,op); }
    void copy (const FabArray<FAB>& src,
               int                  src_comp,
               int                  dest_comp,
               int                  num_comp,
	       const IntVect&       src_nghost,
	       const IntVect&       dst_nghost,
	       const Periodicity&   period = Periodicity::NonPeriodic(),
               CpOp                 op = FabArrayBase::COPY)
        { ParallelCopy(src,src_comp,dest_comp,num_comp,src_nghost,dst_nghost,period,op); }

    //! Copy from src to this.  this and src have the same BoxArray, but different DistributionMapping
    void Redistribute (const FabArray<FAB>& src,
                       int                  src_comp,
                       int                  dest_comp,
                       int                  num_comp,
                       const IntVect&       nghost);

    //
    // In the following copyTo functions, the destination FAB is identical on each process!!
    //

    /**
    * brief Copy the values contained in the intersection of the
    * valid + nghost region of this FabArray with the FAB dest into dest.
    */
    void copyTo (FAB& dest,
		 int  nghost = 0) const;

    /**
    * \brief Copy the values contained in the intersection of the
    * valid + nghost region of this FabArray with the FAB dest and the Box
    * subbox into that subregion of dest.
    */
    void copyTo (FAB&       dest,
                 const Box& subbox,
		 int        nghost = 0) const;

    /**
    * \brief Copy the values contained in the intersection of the
    * num_comp component valid + nghost region of this FabArray, starting at
    * component src_comp, with the FAB dest into dest, starting at
    * component dest_comp in dest.
    */
    void copyTo (FAB& dest,
		 int  src_comp,
		 int  dest_comp,
		 int  num_comp,
		 int  nghost = 0) const;

    /**
    * \brief Copy the values contained in the intersection of the
    * num_comp component valid + nghost region of this FabArray, starting at
    * component src_comp, with the FAB dest and the Box subbox, into
    * dest, starting at component dest_comp in dest.
    */
    void copyTo (FAB&       dest,
		 const Box& subbox,
		 int        src_comp,
		 int        dest_comp,
		 int        num_comp,
		 int        nghost = 0) const;

    //! Shift the boxarray by vector v
    void shift (const IntVect& v);

    bool defined (int i) const noexcept;
    bool defined (const MFIter& mfi) const noexcept;

    /**
    * \brief Copy on intersection within a FabArray.  Data is copied from
    * valid regions to intersecting regions of definition.  The
    * purpose is to fill in the boundary regions of each FAB in
    * the FabArray.  If cross=true, corner cells are not filled.
    * If the length of periodic is provided, periodic boundaries are
    * also filled.  Note that FabArray itself does not contains
    * any periodicity information.
    * FillBoundary expects that its cell-centered version of its BoxArray
    * is non-overlapping.
    */
    void FillBoundary (bool cross = false);

    void FillBoundary (const Periodicity& period, bool cross = false);
    void FillBoundary (const IntVect& nghost, const Periodicity& period, bool cross = false);

    //! Same as FillBoundary(), but only copies ncomp components starting at scomp.
    void FillBoundary (int scomp, int ncomp, bool cross = false);
    void FillBoundary (int scomp, int ncomp, const Periodicity& period, bool cross = false);
    void FillBoundary (int scomp, int ncomp, const IntVect& nghost, const Periodicity& period, bool cross = false);

    void FillBoundary_nowait (bool cross = false);
    void FillBoundary_nowait (const Periodicity& period, bool cross = false);
    void FillBoundary_nowait (int scomp, int ncomp, bool cross = false);
    void FillBoundary_nowait (int scomp, int ncomp, const Periodicity& period, bool cross = false);
    void FillBoundary_nowait (int scomp, int ncomp, const IntVect& nghost, const Periodicity& period, bool cross = false);
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void FillBoundary_finish ();

    void FillBoundary_test ();

    /** \brief Fill cells outside periodic domains with their corresponding cells inside
    * the domain.  Ghost cells are treated the same as valid cells.  The BoxArray
    * is allowed to be overlapping.
    */
    void EnforcePeriodicity (const Periodicity& period);
    void EnforcePeriodicity (int scomp, int ncomp, const Periodicity& period);
    void EnforcePeriodicity (int scomp, int ncomp, const IntVect& nghost,
                             const Periodicity& period);

    // covered   : ghost cells covered by valid cells of this FabArray
    //             (including periodically shifted valid cells)
    // notcovered: ghost cells not covered by valid cells
    //             (including ghost cells outside periodic boundaries)
    // physbnd   : boundary cells outside the domain (excluding periodic boundaries)
    // interior  : interior cells (i.e., valid cells)
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void BuildMask (const Box& phys_domain, const Periodicity& period,
		    value_type covered, value_type notcovered,
		    value_type physbnd, value_type interior);

    // The following are private functions.  But we have to make them public for cuda.

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void FBEP_nowait (int scomp, int ncomp, const IntVect& nghost,
                      const Periodicity& period, bool cross,
                      bool enforce_periodicity_only = false);

    void FB_local_copy_cpu (const FB& TheFB, int scomp, int ncomp);
    void PC_local_cpu (const CPC& thecpc, FabArray<FAB> const& src,
                       int scomp, int dcomp, int ncomp, CpOp op);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    void setVal (value_type x, const CommMetaData& thecmd, int scomp, int ncomp);

    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
    LayoutData<int> RecvLayoutMask (const CommMetaData& thecmd);

#ifdef AMREX_USE_GPU

    void FB_local_copy_gpu (const FB& TheFB, int scomp, int ncomp);
    void PC_local_gpu (const CPC& thecpc, FabArray<FAB> const& src,
                       int scomp, int dcomp, int ncomp, CpOp op);

    void CMD_local_setVal_gpu (value_type x, const CommMetaData& thecmd, int scomp, int ncomp);
    void CMD_remote_setVal_gpu (value_type x, const CommMetaData& thecmd, int scomp, int ncomp);

#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10))

    void FB_local_copy_cuda_graph_1 (const FB& TheFB, int scomp, int ncomp);
    void FB_local_copy_cuda_graph_n (const FB& TheFB, int scomp, int ncomp);

#endif
#endif

#ifdef AMREX_USE_MPI

#ifdef AMREX_USE_GPU
#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10) )

    void FB_pack_send_buffer_cuda_graph (const FB& TheFB, int scomp, int ncomp,
                                         Vector<char*>& send_data,
                                         Vector<std::size_t> const& send_size,
                                         Vector<const CopyComTagsContainer*> const& send_cctc);

    void FB_unpack_recv_buffer_cuda_graph (const FB& TheFB, int dcomp, int ncomp,
                                           Vector<char*> const& recv_data,
                                           Vector<std::size_t> const& recv_size,
                                           Vector<const CopyComTagsContainer*> const& recv_cctc,
                                           bool is_thread_safe);

#endif

    static void pack_send_buffer_gpu (FabArray<FAB> const& src, int scomp, int ncomp,
                                      Vector<char*>& send_data,
                                      Vector<std::size_t> const& send_size,
                                      Vector<const CopyComTagsContainer*> const& send_cctc);

    static void unpack_recv_buffer_gpu (FabArray<FAB>& dst, int dcomp, int ncomp,
                                        Vector<char*> const& recv_data,
                                        Vector<std::size_t> const& recv_size,
                                        Vector<const CopyComTagsContainer*> const& recv_cctc,
                                        CpOp op, bool is_thread_safe);

#endif

    static void pack_send_buffer_cpu (FabArray<FAB> const& src, int scomp, int ncomp,
                                      Vector<char*>& send_data,
                                      Vector<std::size_t> const& send_size,
                                      Vector<const CopyComTagsContainer*> const& send_cctc);

    static void unpack_recv_buffer_cpu (FabArray<FAB>& dst, int dcomp, int ncomp,
                                        Vector<char*> const& recv_data,
                                        Vector<std::size_t> const& recv_size,
                                        Vector<const CopyComTagsContainer*> const& recv_cctc,
                                        CpOp op, bool is_thread_safe);

#endif

protected:

    std::unique_ptr<FabFactory<FAB> > m_factory;
    DataAllocator m_dallocator;

    bool define_function_called = false;

    //
    //! The data.
    std::vector<FAB*> m_fabs_v;

    Vector<std::string> m_tags;

    //! for shared memory
    struct ShMem {
	ShMem () noexcept : alloc(false), n_values(0), n_points(0)
#if defined(BL_USE_MPI3)
		 , win(MPI_WIN_NULL)
#endif
	    { }
	~ShMem () {
#if defined(BL_USE_MPI3)
	    if (win != MPI_WIN_NULL) MPI_Win_free(&win);
#endif
#ifdef BL_USE_TEAM
	    if (alloc) {
		amrex::update_fab_stats(-n_points, -n_values, sizeof(value_type));
            }
#endif
	}
	ShMem (ShMem&& rhs) noexcept
                 : alloc(rhs.alloc), n_values(rhs.n_values), n_points(rhs.n_points)
#if defined(BL_USE_MPI3)
		 , win(rhs.win)
#endif
	{
	    rhs.alloc = false;
#if defined(BL_USE_MPI3)
	    rhs.win = MPI_WIN_NULL;
#endif
	}
	ShMem& operator= (ShMem&& rhs) noexcept {
            if (&rhs != this) {
                alloc = rhs.alloc;
                n_values = rhs.n_values;
                n_points = rhs.n_points;
                rhs.alloc = false;
#if defined(BL_USE_MPI3)
                win = rhs.win;
                rhs.win = MPI_WIN_NULL;
#endif
            }
            return *this;
        }
	ShMem (const ShMem&) = delete;
	ShMem& operator= (const ShMem&) = delete;
	bool  alloc;
	Long  n_values;
	Long  n_points;
#if defined(BL_USE_MPI3)
	MPI_Win win;
#endif
    };
    ShMem shmem;

    bool SharedMemory () const noexcept { return shmem.alloc; }

private:
    typedef typename std::vector<FAB*>::iterator    Iterator;

    void AllocFabs (const FabFactory<FAB>& factory, Arena* ar,
                    const Vector<std::string>& tags);

#ifdef BL_USE_MPI
    //! Prepost nonblocking receives
    void PostRcvs (const MapOfCopyComTagContainers&       m_RcvTags,
                   char*&                                 the_recv_data,
                   Vector<char*>&                         recv_data,
                   Vector<std::size_t>&                   recv_size,
                   Vector<int>&                           recv_from,
                   Vector<MPI_Request>&                   recv_reqs,
                   int                                    ncomp,
                   int                                    SeqNum);

#endif

public:
    //! Data used in non-blocking FillBoundary
    bool fb_cross, fb_epo;
    int fb_scomp, fb_ncomp;
    IntVect fb_nghost;
    Periodicity fb_period;

    //
    char*               fb_the_recv_data = nullptr;
    char*               fb_the_send_data = nullptr;
    Vector<int>         fb_recv_from;
    Vector<char*>       fb_recv_data;
    Vector<std::size_t> fb_recv_size;
    Vector<MPI_Request> fb_recv_reqs;
    Vector<MPI_Status>  fb_recv_stat;
    //
    Vector<char*>       fb_send_data;
    Vector<MPI_Request> fb_send_reqs;
    int                 fb_tag;
};


#include <AMReX_FabArrayCommI.H>

template <class FAB>
bool
FabArray<FAB>::defined (int K) const noexcept
{
    int li = localindex(K);
    if (li >= 0 && li < m_fabs_v.size() && m_fabs_v[li] != 0) {
	return true;
    }
    else {
	return false;
    }
}

template <class FAB>
bool
FabArray<FAB>::defined (const MFIter& mfi) const noexcept
{
    int li = mfi.LocalIndex();
    if (li < static_cast<int>(m_fabs_v.size()) && m_fabs_v[li] != 0) {
	return true;
    }
    else {
	return false;
    }
}

template <class FAB>
FAB*
FabArray<FAB>::fabPtr (const MFIter& mfi) noexcept
{
    BL_ASSERT(mfi.LocalIndex() < indexArray.size());
    BL_ASSERT(DistributionMap() == mfi.DistributionMap());
    int li = mfi.LocalIndex();
    return m_fabs_v[li];
}

template <class FAB>
FAB const*
FabArray<FAB>::fabPtr (const MFIter& mfi) const noexcept
{
    BL_ASSERT(mfi.LocalIndex() < indexArray.size());
    BL_ASSERT(DistributionMap() == mfi.DistributionMap());
    int li = mfi.LocalIndex();
    return m_fabs_v[li];
}

template <class FAB>
FAB*
FabArray<FAB>::fabPtr (int K) noexcept
{
    int li = localindex(K);
    BL_ASSERT(li >=0 && li < indexArray.size());
    return m_fabs_v[li];
}

template <class FAB>
FAB const*
FabArray<FAB>::fabPtr (int K) const noexcept
{
    int li = localindex(K);
    BL_ASSERT(li >=0 && li < indexArray.size());
    return m_fabs_v[li];
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::prefetchToHost (const MFIter& mfi) const noexcept
{
#ifdef AMREX_USE_CUDA
    this->fabPtr(mfi)->prefetchToHost();
#else
    amrex::ignore_unused(mfi);
#endif
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::prefetchToDevice (const MFIter& mfi) const noexcept
{
#ifdef AMREX_USE_CUDA
    this->fabPtr(mfi)->prefetchToDevice();
#else
    amrex::ignore_unused(mfi);
#endif
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::array (const MFIter& mfi) const noexcept
{
    return fabPtr(mfi)->const_array();
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type>
FabArray<FAB>::array (const MFIter& mfi) noexcept
{
    return fabPtr(mfi)->array();
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::array (int K) const noexcept
{
    return fabPtr(K)->const_array();
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type>
FabArray<FAB>::array (int K) noexcept
{
    return fabPtr(K)->array();
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::const_array (const MFIter& mfi) const noexcept
{
    return fabPtr(mfi)->const_array();
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::const_array (int K) const noexcept
{
    return fabPtr(K)->const_array();
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::array (const MFIter& mfi, int start_comp) const noexcept
{
    return fabPtr(mfi)->const_array(start_comp);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type>
FabArray<FAB>::array (const MFIter& mfi, int start_comp) noexcept
{
    return fabPtr(mfi)->array(start_comp);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::array (int K, int start_comp) const noexcept
{
    return fabPtr(K)->const_array(start_comp);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type>
FabArray<FAB>::array (int K, int start_comp) noexcept
{
    return fabPtr(K)->array(start_comp);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::const_array (const MFIter& mfi, int start_comp) const noexcept
{
    return fabPtr(mfi)->const_array(start_comp);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
Array4<typename FabArray<FAB>::value_type const>
FabArray<FAB>::const_array (int K, int start_comp) const noexcept
{
    return fabPtr(K)->const_array(start_comp);
}

template <class FAB>
void
FabArray<FAB>::clear ()
{
    if (define_function_called)
    {
        define_function_called = false;
        clearThisBD();  //!< addThisBD is called in define
    }

    Long nbytes = 0L;
    for (auto x : m_fabs_v) {
        if (x) nbytes += amrex::nBytesOwned(*x);
        m_factory->destroy(x);
    }
    m_fabs_v.clear();
    m_factory.reset();
    m_dallocator.m_arena = nullptr;
    // no need to clear the non-blocking fillboundary stuff

    if (nbytes > 0) {
        for (auto const& t : m_tags) {
            updateMemUsage(t, -nbytes, nullptr);
        }
    }
    m_tags.clear();

    FabArrayBase::clear();
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setVal (value_type val, int nghost)
{
    setVal(val,0,n_comp,IntVect(nghost));
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setVal (value_type val, const IntVect& nghost)
{
    setVal(val,0,n_comp,nghost);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setVal (value_type val, const Box& region, int nghost)
{
    setVal(val,region,0,n_comp,IntVect(nghost));
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setVal (value_type val, const Box& region, const IntVect& nghost)
{
    setVal(val,region,0,n_comp,nghost);
}

template <class FAB>
FabArray<FAB>::FabArray ()
    : shmem()
{
    m_FA_stats.recordBuild();
}

template <class FAB>
FabArray<FAB>::FabArray (const BoxArray&            bxs,
                         const DistributionMapping& dm,
                         int                        nvar,
                         int                        ngrow,
			 const MFInfo&              info,
                         const FabFactory<FAB>&     factory)
    : FabArray<FAB>(bxs,dm,nvar,IntVect(ngrow),info,factory)
{}

template <class FAB>
FabArray<FAB>::FabArray (const BoxArray&            bxs,
                         const DistributionMapping& dm,
                         int                        nvar,
                         const IntVect&             ngrow,
			 const MFInfo&              info,
                         const FabFactory<FAB>&     factory)
    : m_factory(factory.clone()),
      shmem()
{
    m_FA_stats.recordBuild();
    define(bxs,dm,nvar,ngrow,info,*m_factory);
}

template <class FAB>
FabArray<FAB>::FabArray (const FabArray<FAB>& rhs, MakeType maketype, int scomp, int ncomp)
    : m_factory(rhs.Factory().clone()),
      shmem()
{
    m_FA_stats.recordBuild();
    define(rhs.boxArray(), rhs.DistributionMap(), ncomp, rhs.nGrowVect(),
           MFInfo().SetAlloc(false), *m_factory);

    if (maketype == amrex::make_alias)
    {
        for (int i = 0, n = indexArray.size(); i < n; ++i) {
            auto const& rhsfab = *(rhs.m_fabs_v[i]);
            m_fabs_v.push_back(m_factory->create_alias(rhsfab, scomp, ncomp));
        }
    }
    else
    {
        amrex::Abort("FabArray: unknown MakeType");
    }
}

template <class FAB>
FabArray<FAB>::FabArray (FabArray<FAB>&& rhs) noexcept
    : FabArrayBase (std::move(rhs))
    , m_factory    (std::move(rhs.m_factory))
    , m_dallocator (std::move(rhs.m_dallocator))
    , define_function_called(rhs.define_function_called)
    , m_fabs_v     (std::move(rhs.m_fabs_v))
    , m_tags       (std::move(rhs.m_tags))
    , shmem        (std::move(rhs.shmem))
    // no need to worry about the data used in non-blocking FillBoundary.
{
    m_FA_stats.recordBuild();
    rhs.define_function_called = false; // the responsibility of clear BD has been transferred.
    rhs.m_fabs_v.clear(); // clear the data pointers so that rhs.clear does delete them.
    rhs.clear();
}

template <class FAB>
FabArray<FAB>&
FabArray<FAB>::operator= (FabArray<FAB>&& rhs) noexcept
{
    if (&rhs != this)
    {
        clear();

        FabArrayBase::operator=(std::move(rhs));
        m_factory = std::move(rhs.m_factory);
        m_dallocator = std::move(rhs.m_dallocator);
        define_function_called = rhs.define_function_called;
        std::swap(m_fabs_v, rhs.m_fabs_v);
        std::swap(m_tags, rhs.m_tags);
        shmem = std::move(rhs.shmem);

        rhs.define_function_called = false;
        rhs.m_fabs_v.clear();
        rhs.m_tags.clear();
        rhs.clear();
    }
    return *this;
}

template <class FAB>
FabArray<FAB>::~FabArray ()
{
    m_FA_stats.recordDelete();
    clear();
}

template <class FAB>
bool
FabArray<FAB>::ok () const
{
    if (!define_function_called) return false;

    int isok = 1;

    for (MFIter fai(*this); fai.isValid() && isok; ++fai)
    {
        if (defined(fai))
        {
            if (get(fai).box() != fabbox(fai.index()))
            {
                isok = 0;
            }
        }
        else
        {
            isok = 0;
        }
    }

    ParallelAllReduce::Min(isok, ParallelContext::CommunicatorSub());

    return isok == 1;
}

template <class FAB>
void
FabArray<FAB>::define (const BoxArray&            bxs,
		       const DistributionMapping& dm,
		       int                        nvar,
		       int                        ngrow,
		       const MFInfo&              info,
                       const FabFactory<FAB>&     a_factory)
{
    define(bxs,dm,nvar,IntVect(ngrow),info,a_factory);
}

template <class FAB>
void
FabArray<FAB>::define (const BoxArray&            bxs,
		       const DistributionMapping& dm,
		       int                        nvar,
		       const IntVect&             ngrow,
		       const MFInfo&              info,
                       const FabFactory<FAB>&     a_factory)
{
    std::unique_ptr<FabFactory<FAB> > factory(a_factory.clone());

    clear();

    m_factory = std::move(factory);
    m_dallocator.m_arena = info.arena;

    define_function_called = true;

    BL_ASSERT(ngrow.allGE(IntVect::TheZeroVector()));
    BL_ASSERT(boxarray.size() == 0);
    FabArrayBase::define(bxs, dm, nvar, ngrow);

    addThisBD();

    if(info.alloc) {
        AllocFabs(*m_factory, info.arena, info.tags);
        Gpu::synchronize();
#ifdef BL_USE_TEAM
        ParallelDescriptor::MyTeam().MemoryBarrier();
#endif
    }
}

template <class FAB>
void
FabArray<FAB>::AllocFabs (const FabFactory<FAB>& factory, Arena* ar,
                          const Vector<std::string>& tags)
{
    const int n = indexArray.size();
    const int nworkers = ParallelDescriptor::TeamSize();
    shmem.alloc = (nworkers > 1);

    bool alloc = !shmem.alloc;

    FabInfo fab_info;
    fab_info.SetAlloc(alloc).SetShared(shmem.alloc).SetArena(ar);

    m_fabs_v.reserve(n);

    Long nbytes = 0L;
    for (int i = 0; i < n; ++i)
    {
        int K = indexArray[i];
        const Box& tmpbox = fabbox(K);
        m_fabs_v.push_back(factory.create(tmpbox, n_comp, fab_info, K));
        nbytes += amrex::nBytesOwned(*m_fabs_v.back());
    }

    m_tags.clear();
    m_tags.emplace_back("All");
    for (auto const& t : m_region_tag) {
        m_tags.push_back(t);
    }
    for (auto const& t : tags) {
        m_tags.push_back(t);
    }
    for (auto const& t: m_tags) {
        updateMemUsage(t, nbytes, ar);
    }

#ifdef BL_USE_TEAM
    if (shmem.alloc)
    {
        const int teamlead = ParallelDescriptor::MyTeamLead();

        shmem.n_values = 0;
        shmem.n_points = 0;
        Vector<Long> offset(n,0);
        Vector<Long> nextoffset(nworkers,-1);
        for (int i = 0; i < n; ++i) {
            int K = indexArray[i];
            int owner = distributionMap[K] - teamlead;
            Long s = m_fabs_v[i]->size();
            if (ownership[i]) {
                shmem.n_values += s;
                shmem.n_points += m_fabs_v[i]->numPts();
            }
            if (nextoffset[owner] < 0) {
                offset[i] = 0;
                nextoffset[owner] = s;
            } else {
                offset[i] = nextoffset[owner];
                nextoffset[owner] += s;
            }
        }

        size_t bytes = shmem.n_values*sizeof(value_type);

        value_type *mfp;
        Vector<value_type*> dps;

#if defined (BL_USE_MPI3)

        static MPI_Info info = MPI_INFO_NULL;
        if (info == MPI_INFO_NULL) {
            MPI_Info_create(&info);
            MPI_Info_set(info, "alloc_shared_noncontig", "true");
        }

        const MPI_Comm& team_comm = ParallelDescriptor::MyTeam().get();

        BL_MPI_REQUIRE( MPI_Win_allocate_shared(bytes, sizeof(value_type),
                                                info, team_comm, &mfp, &shmem.win) );

        for (int w = 0; w < nworkers; ++w) {
            MPI_Aint sz;
            int disp;
            value_type *dptr = 0;
            BL_MPI_REQUIRE( MPI_Win_shared_query(shmem.win, w, &sz, &disp, &dptr) );
            // BL_ASSERT(disp == sizeof(value_type));
            dps.push_back(dptr);
        }

#else

        amrex::Abort("BaseFab::define: to allocate shared memory, USE_MPI3 must be true");

#endif

        for (int i = 0; i < n; ++i) {
            int K = indexArray[i];
            int owner = distributionMap[K] - teamlead;
            value_type *p = dps[owner] + offset[i];
            m_fabs_v[i]->setPtr(p, m_fabs_v[i]->size());
        }

        for (Long i = 0; i < shmem.n_values; i++, mfp++) {
            new (mfp) value_type;
        }

        amrex::update_fab_stats(shmem.n_points, shmem.n_values, sizeof(value_type));
    }
#endif
}

template <class FAB>
void
FabArray<FAB>::setFab (int  boxno,
                       FAB* elem)
{
    //
    // Must check it is of the proper size.
    //
    if (n_comp == 0) {
        n_comp = elem->nComp();
    }

    BL_ASSERT(n_comp == elem->nComp());
    BL_ASSERT(boxarray.size() > 0);
    BL_ASSERT(elem->box() == fabbox(boxno));
    BL_ASSERT(!this->defined(boxno));
    BL_ASSERT(distributionMap[boxno] == ParallelDescriptor::MyProc());

    if (m_fabs_v.size() == 0) {
      m_fabs_v.resize(indexArray.size(),nullptr);
    }

    const int li = localindex(boxno);
    m_fabs_v[li] = elem;
}

template <class FAB>
void
FabArray<FAB>::setFab (const MFIter& mfi,
                       FAB* elem, bool assertion)
{
    if (assertion)
    {
        //
        // Must check it is of the proper size.
        //
        if (n_comp == 0)
            n_comp = elem->nComp();

        BL_ASSERT(n_comp == elem->nComp());
        BL_ASSERT(boxarray.size() > 0);
        BL_ASSERT(elem->box() == fabbox(mfi.index()));
        BL_ASSERT(!this->defined(mfi));
        BL_ASSERT(distributionMap[mfi.index()] == ParallelDescriptor::MyProc());
    }

    if (m_fabs_v.size() == 0) {
      m_fabs_v.resize(indexArray.size(),nullptr);
    }

    const int li = mfi.LocalIndex();
    m_fabs_v[li] = elem;
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setBndry (value_type val)
{
    setBndry(val, 0, n_comp);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setBndry (value_type val,
                         int        strt_comp,
                         int        ncomp)
{
    if (n_grow.max() > 0)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter fai(*this); fai.isValid(); ++fai)
        {
            get(fai).template setComplement<RunOn::Device>(val, fai.validbox(), strt_comp, ncomp);
        }
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setDomainBndry (value_type val, const Geometry& geom)
{
    setDomainBndry(val, 0, n_comp, geom);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setDomainBndry (value_type val,
                               int        strt_comp,
                               int        ncomp,
                               const Geometry& geom)
{
    Box domain_box = amrex::convert(geom.Domain(), boxArray().ixType());
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) {
            int n = domain_box.length(idim);
            domain_box.grow(idim, n);
        }
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        const Box& gbx = fai.fabbox();
        if (not domain_box.contains(gbx))
        {
            get(fai).template setComplement<RunOn::Device>(val, domain_box, strt_comp, ncomp);
        }
    }
}

//
// Copies to FABs, note that destination is first arg.
//

template <class FAB>
void
FabArray<FAB>::copyTo (FAB& dest,
		       int  nghost) const
{
    copyTo(dest, dest.box(), 0, 0, dest.nComp(), nghost);
}

template <class FAB>
void
FabArray<FAB>::copyTo (FAB&       dest,
		       const Box& subbox,
		       int        nghost) const
{
    copyTo(dest, subbox, 0, 0, dest.nComp(), nghost);
}

template <class FAB>
void
FabArray<FAB>::copyTo (FAB& dest,
		       int  scomp,
		       int  dcomp,
		       int  ncomp,
		       int  nghost) const
{
    copyTo(dest, dest.box(), scomp, dcomp, ncomp, nghost);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setVal (value_type val)
{
    setVal(val,0,n_comp,n_grow);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::operator= (value_type val)
{
    setVal(val);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setVal (value_type val,
                       int        comp,
                       int        ncomp,
                       int        nghost)
{
    setVal(val,comp,ncomp,IntVect(nghost));
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::setVal (value_type val,
                       int        comp,
                       int        ncomp,
                       const IntVect& nghost)
{
    BL_ASSERT(nghost.allGE(IntVect::TheZeroVector()) && nghost.allLE(n_grow));
    BL_ASSERT(comp+ncomp <= n_comp);

    BL_PROFILE("FabArray::setVal()");

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter fai(*this,TilingIfNotGPU()); fai.isValid(); ++fai)
    {
	const Box& bx = fai.growntilebox(nghost);
        auto fab = this->array(fai);
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            fab(i,j,k,n+comp) = val;
        });
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setVal (value_type val,
                       const Box& region,
                       int        comp,
                       int        ncomp,
                       int        nghost)
{
    setVal(val,region,comp,ncomp,IntVect(nghost));
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::setVal (value_type val,
                       const Box& region,
                       int        comp,
                       int        ncomp,
                       const IntVect& nghost)
{
    BL_ASSERT(nghost.allGE(IntVect::TheZeroVector()) && nghost.allLE(n_grow));
    BL_ASSERT(comp+ncomp <= n_comp);

    BL_PROFILE("FabArray::setVal(val,region,comp,ncomp,nghost)");

#ifdef _OPENMP
    AMREX_ALWAYS_ASSERT(!omp_in_parallel());
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter fai(*this,TilingIfNotGPU()); fai.isValid(); ++fai)
    {
        Box b = fai.growntilebox(nghost) & region;

        if (b.ok()) {
            auto fab = this->array(fai);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( b, ncomp, i, j, k, n,
            {
                fab(i,j,k,n+comp) = val;
            });
        }
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::abs (int comp, int ncomp, int nghost)
{
    abs(comp, ncomp, IntVect(nghost));
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::abs (int comp, int ncomp, const IntVect& nghost)
{
    BL_ASSERT(nghost.allGE(IntVect::TheZeroVector()) && nghost.allLE(n_grow));
    BL_ASSERT(comp+ncomp <= n_comp);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        auto fab = this->array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            fab(i,j,k,n+comp) = amrex::Math::abs(fab(i,j,k,n+comp));
        });
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::plus (value_type val, int comp, int num_comp, int nghost)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        auto fab = this->array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, num_comp, i, j, k, n,
        {
            fab(i,j,k,n+comp) += val;
        });
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::plus (value_type val, const Box& region, int comp, int num_comp, int nghost)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost) & region;
        if (bx.ok()) {
            auto fab = this->array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, num_comp, i, j, k, n,
            {
                fab(i,j,k,n+comp) += val;
            });
        }
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::mult (value_type val, int comp, int num_comp, int nghost)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        auto fab = this->array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, num_comp, i, j, k, n,
        {
            fab(i,j,k,n+comp) *= val;
        });
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::mult (value_type val, const Box& region, int comp, int num_comp, int nghost)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost) & region;
        if (bx.ok()) {
            auto fab = this->array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, num_comp, i, j, k, n,
            {
                fab(i,j,k,n+comp) *= val;
            });
        }
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::invert (value_type numerator, int comp, int num_comp, int nghost)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        auto fab = this->array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, num_comp, i, j, k, n,
        {
            fab(i,j,k,n+comp) = numerator / fab(i,j,k,n+comp);
        });
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::invert (value_type numerator, const Box& region, int comp, int num_comp, int nghost)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost) & region;
        if (bx.ok()) {
            auto fab = this->array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, num_comp, i, j, k, n,
            {
                fab(i,j,k,n+comp) = numerator / fab(i,j,k,n+comp);
            });
        }
    }
}

template <class FAB>
void
FabArray<FAB>::shift (const IntVect& v)
{
    clearThisBD();  // The new boxarry will have a different ID.
    for(int id(0); id < AMREX_SPACEDIM; ++id)
    {
      boxarray.shift(id, v[id]);
    }
    addThisBD();
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        get(fai).shift(v);
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundary (bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    if ( n_grow.max() > 0 ) {
	FillBoundary_nowait(0, nComp(), n_grow, Periodicity::NonPeriodic(), cross);
	FillBoundary_finish();
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundary (const Periodicity& period, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    if ( n_grow.max() > 0 ) {
	FillBoundary_nowait(0, nComp(), n_grow, period, cross);
	FillBoundary_finish();
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundary (const IntVect& nghost, const Periodicity& period, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nghost.allLE(nGrowVect()),
                                     "FillBoundary: asked to fill more ghost cells than we have");
    if ( nghost.max() > 0 ) {
	FillBoundary_nowait(0, nComp(), nghost, period, cross);
	FillBoundary_finish();
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundary (int scomp, int ncomp, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    if ( n_grow.max() > 0 ) {
	FillBoundary_nowait(scomp, ncomp, n_grow, Periodicity::NonPeriodic(), cross);
	FillBoundary_finish();
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundary (int scomp, int ncomp, const Periodicity& period, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    if ( n_grow.max() > 0 ) {
	FillBoundary_nowait(scomp, ncomp, n_grow, period, cross);
	FillBoundary_finish();
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundary (int scomp, int ncomp, const IntVect& nghost,
                             const Periodicity& period, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nghost.allLE(nGrowVect()),
                                     "FillBoundary: asked to fill more ghost cells than we have");
    if ( nghost.max() > 0 ) {
	FillBoundary_nowait(scomp, ncomp, nghost, period, cross);
	FillBoundary_finish();
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundary_nowait (bool cross)
{
    FillBoundary_nowait(0, nComp(), nGrowVect(), Periodicity::NonPeriodic(), cross);
}

template <class FAB>
void
FabArray<FAB>::FillBoundary_nowait (const Periodicity& period, bool cross)
{
    FillBoundary_nowait(0, nComp(), nGrowVect(), period, cross);
}

template <class FAB>
void
FabArray<FAB>::FillBoundary_nowait (int scomp, int ncomp, bool cross)
{
    FillBoundary_nowait(scomp, ncomp, nGrowVect(), Periodicity::NonPeriodic(), cross);
}

template <class FAB>
void
FabArray<FAB>::EnforcePeriodicity (const Periodicity& period)
{
    BL_PROFILE("FabArray::EnforcePeriodicity");
    if (period.isAnyPeriodic()) {
	FBEP_nowait(0, nComp(), nGrowVect(), period, false, true);
	FillBoundary_finish(); // unsafe unless isAnyPeriodic()
    }
}

template <class FAB>
void
FabArray<FAB>::EnforcePeriodicity (int scomp, int ncomp, const Periodicity& period)
{
    BL_PROFILE("FabArray::EnforcePeriodicity");
    if (period.isAnyPeriodic()) {
	FBEP_nowait(scomp, ncomp, nGrowVect(), period, false, true);
	FillBoundary_finish(); // unsafe unless isAnyPeriodic()
    }
}

template <class FAB>
void
FabArray<FAB>::EnforcePeriodicity (int scomp, int ncomp, const IntVect& nghost,
                                   const Periodicity& period)
{
    BL_PROFILE("FabArray::EnforcePeriodicity");
    if (period.isAnyPeriodic()) {
	FBEP_nowait(scomp, ncomp, nghost, period, false, true);
	FillBoundary_finish(); // unsafe unless isAnyPeriodic()
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundary_nowait (int scomp, int ncomp, const Periodicity& period, bool cross)
{
    BL_PROFILE("FillBoundary_nowait()");
    FBEP_nowait(scomp, ncomp, nGrowVect(), period, cross);
}

template <class FAB>
void
FabArray<FAB>::FillBoundary_nowait (int scomp, int ncomp, const IntVect& nghost,
                                    const Periodicity& period, bool cross)
{
    BL_PROFILE("FillBoundary_nowait()");
    FBEP_nowait(scomp, ncomp, nghost, period, cross);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type Z>
void
FabArray<FAB>::BuildMask (const Box& phys_domain, const Periodicity& period,
			  value_type covered, value_type notcovered,
			  value_type physbnd, value_type interior)
{
    int ncomp = this->nComp();
    const IntVect& ngrow = this->nGrowVect();

    Box domain = amrex::convert(phys_domain, boxArray().ixType());
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
	if (period.isPeriodic(i)) {
	    domain.grow(i, ngrow[i]);
	}
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& fab = this->array(mfi);
        Box const& fbx = mfi.growntilebox();
        Box const& gbx = fbx & domain;
        Box const& vbx = mfi.validbox();
        AMREX_HOST_DEVICE_FOR_4D(fbx, ncomp, i, j, k, n,
        {
            IntVect iv(AMREX_D_DECL(i,j,k));
            if (vbx.contains(iv)) {
                fab(i,j,k,n) = interior;
            } else if (gbx.contains(iv)) {
                fab(i,j,k,n) = notcovered;
            } else {
                fab(i,j,k,n) = physbnd;
            }
        });
    }

    const FabArrayBase::FB& TheFB = this->getFB(ngrow,period);
    setVal(covered, TheFB, 0, ncomp);
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
void
FabArray<FAB>::setVal (value_type val, const CommMetaData& thecmd, int scomp, int ncomp)
{
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        CMD_local_setVal_gpu(val, thecmd, scomp, ncomp);
        CMD_remote_setVal_gpu(val, thecmd, scomp, ncomp);
    }
    else
#endif
    {
        const CopyComTagsContainer&      LocTags = *(thecmd.m_LocTags);
        const MapOfCopyComTagContainers& RcvTags = *(thecmd.m_RcvTags);
        int N_locs = LocTags.size();
#ifdef _OPENMP
#pragma omp parallel for if (thecmd.m_threadsafe_loc)
#endif
        for (int i = 0; i < N_locs; ++i) {
            const CopyComTag& tag = LocTags[i];
            (*this)[tag.dstIndex].template setVal<RunOn::Host>(val, tag.dbox, scomp, ncomp);
        }

        for (auto it = RcvTags.begin(); it != RcvTags.end(); ++it) {
            int N = it->second.size();
#ifdef _OPENMP
#pragma omp parallel for if (thecmd.m_threadsafe_rcv)
#endif
            for (int i = 0; i < N; ++i) {
                const CopyComTag& tag = it->second[i];
                (*this)[tag.dstIndex].template setVal<RunOn::Host>(val, tag.dbox, scomp, ncomp);
            }
        }
    }
}

template <class FAB>
template <class F, typename std::enable_if<IsBaseFab<F>::value,int>::type>
LayoutData<int>
FabArray<FAB>::RecvLayoutMask (const CommMetaData& thecmd)
{
    LayoutData<int> r(this->boxArray(), this->DistributionMap());
#ifdef _OPENMP
#pragma omp parallel if (thecmd.m_threadsafe_rcv)
#endif
    for (MFIter mfi(r); mfi.isValid(); ++mfi) {
        r[mfi] = 0;
    }

    const CopyComTagsContainer&      LocTags = *(thecmd.m_LocTags);
    const MapOfCopyComTagContainers& RcvTags = *(thecmd.m_RcvTags);

    int N_locs = LocTags.size();
    for (int i = 0; i < N_locs; ++i) {
        const CopyComTag& tag = LocTags[i];
        r[tag.dstIndex] = 1;
    }

    for (auto it = RcvTags.begin(); it != RcvTags.end(); ++it) {
        int N = it->second.size();
        for (int i = 0; i < N; ++i) {
            const CopyComTag& tag = it->second[i];
            r[tag.dstIndex] = 1;
        }
    }
    return r;
}

}

#endif /*BL_FABARRAY_H*/
