
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

#ifdef BL_USE_UPCXX
#include <AMReX_BLPgas.H>
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
#include <iostream>
#include <AMReX_FabArrayBase.H>
#include <AMReX_MFIter.H>
#include <AMReX_MakeType.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_parstream.H>
#include <AMReX_SPMD.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif

namespace amrex {

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
    MFInfo& SetAlloc(bool a) { alloc = a; return *this; }
};

    template <class T>
    class MFGraph;

template <class FAB>
class FabArray
    :
    public FabArrayBase
{
public:
    friend class Action;
    template <class T>
    friend class MFGraph;

    struct FABType {
        typedef FAB value_type;
    };

    // if FAB is a BaseFab or its child, value_type = FAB::value_type
    // else                              value_type = FAB;
    using value_type = typename std::conditional<IsBaseFab<FAB>::value, FAB, FABType>::type::value_type;

    //
    // Constructs an empty FabArray<FAB>.
    //
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

    FabArray (const FabArray<FAB>& rhs, MakeType maketype, int scomp, int ncomp);

    template <class = typename std::enable_if<IsBaseFab<FAB>::value> >
    FabArray (const BoxArray&            bxs,
              const DistributionMapping& dm,
              int                        nvar,
              int                        ngrow,
              const Vector<value_type*>&  pval);

    //! The destructor -- deletes all FABs in the array.
    virtual ~FabArray ();

    FabArray (FabArray<FAB>&& rhs) noexcept;

    FabArray (const FabArray<FAB>& rhs) = delete;
    FabArray<FAB>& operator= (const FabArray<FAB>& rhs) = delete;
    FabArray<FAB>& operator= (FabArray<FAB>&& rhs) = delete;

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

    const FabFactory<FAB>& Factory () const { return *m_factory; }

    Box getDomain () const { return m_factory->getDomain(); }

    /**
    * \brief Return true if the FabArray is well-defined.  That is,
    * if FABs are allocated for each Box in the BoxArray and the
    * sizes of the FABs and the number of components are consistent
    * with the definition of the FabArray.
    */
    bool ok () const;

    //! Return a constant reference to the FAB associated with mfi.
    const FAB& operator[] (const MFIter& mfi) const;

    //! Return a constant reference to the FAB associated with mfi.
    const FAB& get (const MFIter& mfi) const { return this->operator[](mfi); }

    //! Returns a reference to the FAB associated mfi.
    FAB& operator[] (const MFIter& mfi);

    //! Returns a reference to the FAB associated mfi.
    FAB& get (const MFIter& mfi) { return this->operator[](mfi); }

    //! Return a constant reference to the FAB associated with the Kth element.
    const FAB& operator[] (int K) const;

    //! Return a constant reference to the FAB associated with the Kth element.
    const FAB& get (int K) const { return this->operator[](K); }

    //! Return a reference to the FAB associated with the Kth element.
    FAB& operator[] (int K);

    //! Return a reference to the FAB associated with the Kth element.
    FAB& get (int K)  { return this->operator[](K); }

    //! Explicitly set the Kth FAB in the FabArray to point to elem.
    void setFab (int K, FAB* elem);

    //! Explicitly set the FAB associated with mfi in the FabArray to point to elem.
    void setFab (const MFIter&mfi, FAB* elem, bool assertion=true);

    //! Releases FAB memory in the FabArray.
    void clear ();

    //! Set all components in the entire region of each FAB to val.
    template <class = typename std::enable_if<IsBaseFab<FAB>::value> >
    void setVal (value_type val);

    //! Set all components in the entire region of each FAB to val.
    template <class = typename std::enable_if<IsBaseFab<FAB>::value> >
    void operator= (const value_type& val);

    /**
    * \brief Set the value of num_comp components in the valid region of
    * each FAB in the FabArray, starting at component comp to val.
    * Also set the value of nghost boundary cells.
    */
    template <class = typename std::enable_if<IsBaseFab<FAB>::value> >
    void setVal (value_type val,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);

    /**
    * \brief Set the value of num_comp components in the valid region of
    * each FAB in the FabArray, starting at component comp, as well
    * as nghost boundary cells, to val, provided they also intersect
    * with the Box region.
    */
    template <class = typename std::enable_if<IsBaseFab<FAB>::value> >
    void setVal (value_type val,
                 const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);
    /**
    * \brief Set all components in the valid region of each FAB in the
    * FabArray to val, including nghost boundary cells.
    */
    template <class = typename std::enable_if<IsBaseFab<FAB>::value> >
    void setVal (value_type val,
                 int        nghost);

    /**
    * \brief Set all components in the valid region of each FAB in the
    * FabArray to val, including nghost boundary cells, that also
    * intersect the Box region.
    */
    template <class = typename std::enable_if<IsBaseFab<FAB>::value> >
    void setVal (value_type val,
                 const Box& region,
                 int        nghost);

    //! Set all values in the boundary region to val.
    template <class = typename std::enable_if<IsBaseFab<FAB>::value> >
    void setBndry (value_type val);

    //! Set ncomp values in the boundary region, starting at start_comp to val.
    template <class = typename std::enable_if<IsBaseFab<FAB>::value> >
    void setBndry (value_type val,
                   int        strt_comp,
                   int        ncomp);
 
   //! Set all values outside the Geometry domain to val.
    template <class = typename std::enable_if<IsBaseFab<FAB>::value> >
    void setDomainBndry (value_type val, const Geometry& goem);

    //! Set ncomp values outside the Geometry domain to val, starting at start_comp.
    template <class = typename std::enable_if<IsBaseFab<FAB>::value> >
    void setDomainBndry (value_type val, int strt_comp, int ncomp, const Geometry& goem);
    /**
    * \brief This function copies data from fa to this FabArray.  Each FAB
    * in fa is intersected with all FABs in this FabArray and a copy
    * is performed on the region of intersection.  The intersection
    * is restricted to the valid regions.
    */
    void ParallelCopy (const FabArray<FAB>& fa,
                       const Periodicity&   period = Periodicity::NonPeriodic(),
                       CpOp                 op = FabArrayBase::COPY);
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
    void ParallelCopy (const FabArray<FAB>& src,
                       int                  src_comp,
                       int                  dest_comp,
                       int                  num_comp,
                       const Periodicity&   period = Periodicity::NonPeriodic(),
                       CpOp                 op = FabArrayBase::COPY);
    void copy (const FabArray<FAB>& src,
               int                  src_comp,
               int                  dest_comp,
               int                  num_comp,
	       const Periodicity&   period = Periodicity::NonPeriodic(),
               CpOp                 op = FabArrayBase::COPY)
        { ParallelCopy(src,src_comp,dest_comp,num_comp, period, op); }

    //! Similar to the above function, except that source and destination are grown by src_nghost and dst_nghost, respectively 
    void ParallelCopy (const FabArray<FAB>& src,
                       int                  src_comp,
                       int                  dest_comp,
                       int                  num_comp,
                       int                  src_nghost,
                       int                  dst_nghost,
                       const Periodicity&   period = Periodicity::NonPeriodic(),
                       CpOp                 op = FabArrayBase::COPY);
    void copy (const FabArray<FAB>& src,
               int                  src_comp,
               int                  dest_comp,
               int                  num_comp,
	       int                  src_nghost,
	       int                  dst_nghost,
	       const Periodicity&   period = Periodicity::NonPeriodic(),
               CpOp                 op = FabArrayBase::COPY)
        { ParallelCopy(src,src_comp,dest_comp,num_comp,src_nghost,dst_nghost,period,op); }

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

    bool defined (int i) const;
    bool defined (const MFIter& mfi) const;

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

    //! Same as FillBoundary(), but only copies ncomp components starting at scomp.
    void FillBoundary (int scomp, int ncomp, bool cross = false);
    void FillBoundary (int scomp, int ncomp, const Periodicity& period, bool cross = false);

    void FillBoundary_nowait (bool cross = false);
    void FillBoundary_nowait (const Periodicity& period, bool cross = false);
    void FillBoundary_nowait (int scomp, int ncomp, bool cross = false);
    void FillBoundary_nowait (int scomp, int ncomp, const Periodicity& period, bool cross = false);
    void FillBoundary_finish ();

    /** \brief Fill cells outside periodic domains with their corresponding cells inside
    * the domain.  Ghost cells are treated the same as valid cells.  The BoxArray
    * is allowed to be overlapping.
    */
    void EnforcePeriodicity (const Periodicity& period);
    void EnforcePeriodicity (int scomp, int ncomp, const Periodicity& period);

    // covered   : ghost cells covered by valid cells of this FabArray
    //             (including periodically shifted valid cells)
    // notcovered: ghost cells not covered by valid cells
    //             (including ghost cells outside periodic boundaries)
    // physbnd   : boundary cells outside the domain (excluding periodic boundaries)
    // interior  : interior cells (i.e., valid cells)
    template <class = typename std::enable_if<IsBaseFab<FAB>::value> >
    void BuildMask (const Box& phys_domain, const Periodicity& period,
		    value_type covered, value_type notcovered, 
		    value_type physbnd, value_type interior);

    //! Move FABs in this FabArray to different MPI ranks.
    struct FABMoves {
      int distMapIndex, fromRank, toRank, seqNum;
    };

    int MoveFabs (const Vector<int> &newDistMapArray);
    //static void MoveAllFabs (const Vector<int> &newDistMapArray);
    static void MoveAllFabs (const DistributionMapping &newDistMap);
    static void MoveAllFabs (const DistributionMapping &oldDistMap, const DistributionMapping &newDistMap);
    static void LockAllFAPointers ();
    static void CheckFAPointers (bool abortOnError = true);
    static void PrintFAPointers ();
    virtual void AddProcsToComp (int ioProcNumSCS, int ioProcNumAll,
                                 int scsMyId, MPI_Comm scsComm);

    static void copyInter (FabArray<FAB> *src, FabArray<FAB> *dest,
                           int        src_comp,
                           int        dest_comp,
                           int        num_comp,
                           int        src_nghost,
                           int        dst_nghost,
                           const MPI_Comm &commSrc,
                           const MPI_Comm &commDest,
                           const MPI_Comm &commInter,
                           const MPI_Comm &commBoth,
                           bool       isSrc,
                           CpOp       op = FabArrayBase::COPY);


protected:

    std::unique_ptr<FabFactory<FAB> > m_factory;

    bool define_function_called = false;
    
    //
    // The data.
    //
    std::vector<FAB*> m_fabs_v;

    // ---- currently active fabarrays:  <ngrids to find distmap, <aFAPId, FabArray *> >
    static std::map<int, std::map<int, FabArray<FAB> *> > allocatedFAPointers;

    // ---- set of ids of fabarrays created without allocating memory
    static std::set<int> noallocFAPIds;

    // for shared memory
    struct ShMem {
	ShMem () : alloc(false), n_values(0), n_points(0)
#ifdef BL_USE_UPCXX
		 , p(nullptr)
#elif defined(BL_USE_MPI3)
		 , win(MPI_WIN_NULL)
#endif
	    { }
	~ShMem () {
#ifdef BL_USE_UPCXX
	    if (p) BLPgas::free(p);
#elif defined(BL_USE_MPI3)
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
#ifdef BL_USE_UPCXX
		 , p(rhs.p)
#elif defined(BL_USE_MPI3)
		 , win(rhs.win)
#endif
	{
	    rhs.alloc = false;
#ifdef BL_USE_UPCXX
	    rhs.p = nullptr;
#elif defined(BL_USE_MPI3)
	    rhs.win = MPI_WIN_NULL;
#endif
	}
	ShMem (const ShMem&) = delete;
	ShMem& operator= (const ShMem&) = delete;
	ShMem& operator= (ShMem&&) = delete;
	bool  alloc;
	long  n_values;
	long  n_points;
#ifdef BL_USE_UPCXX
	void *p;
#elif defined(BL_USE_MPI3)
	MPI_Win win;
#endif
    };
    ShMem shmem;

    bool SharedMemory () const { return shmem.alloc; }

private:
    typedef typename std::vector<FAB*>::iterator    Iterator;

    void AllocFabs (const FabFactory<FAB>& factory);

    void FBEP_nowait (int scomp, int ncomp, const Periodicity& period, bool cross,
		      bool enforce_periodicity_only = false);

#ifdef BL_USE_MPI
    //! Prepost nonblocking receives
    void PostRcvs (const MapOfCopyComTagContainers&       m_RcvVols,
                   const MapOfCopyComTagContainers&       m_RcvTags,
                   Vector<char*>&                          recv_data,
                   Vector<int>&                            recv_size,
                   Vector<int>&                            recv_from,
                   Vector<MPI_Request>&                    recv_reqs,
                   int                                    icomp,
                   int                                    ncomp,
                   int                                    SeqNum,
                   int                                    preSeqNum,
                   MPI_Comm   comm = ParallelDescriptor::Communicator());
#endif
    
#ifdef BL_USE_MPI3
    void PostRcvs_MPI_Onesided (const MapOfCopyComTagContainers&  m_RcvVols,
                                char*&                            the_recv_data,
                                Vector<char*>&                     recv_data,
                                Vector<int>&                       recv_size,
                                Vector<int>&                       recv_from,
                                Vector<MPI_Request>&               recv_reqs,
                                Vector<MPI_Aint>          &        recv_disp,
                                int                               icomp,
                                int                               ncomp,
                                int                               SeqNum,
                                MPI_Win &                         win);
#endif
  
#ifdef BL_USE_UPCXX
    void PostRcvs_PGAS (const MapOfCopyComTagContainers&  m_RcvVols,
                        char*&                            the_recv_data,
                        Vector<char*>&                     recv_data,
                        Vector<int>&                       recv_size,
                        Vector<int>&                       recv_from,
                        int                               icomp,
                        int                               ncomp,
                        int                               SeqNum,
                        upcxx::event*                     recv_event);
#endif

public:
    // Data used in non-blocking FillBoundary
    bool fb_cross, fb_epo;
    int fb_scomp, fb_ncomp;
    Periodicity fb_period;

    //
    char*              fb_the_recv_data;
    Vector<int>         fb_recv_from;
    Vector<char*>       fb_recv_data;
    Vector<int>         fb_recv_size;
    Vector<MPI_Request> fb_recv_reqs;
#ifdef BL_USE_MPI3
    Vector<MPI_Aint>    fb_recv_disp;
#endif
    //
    Vector<char*>       fb_send_data;
    Vector<MPI_Request> fb_send_reqs;
    int                fb_tag;
};

#ifdef BL_USE_MPI

template <class FAB>
void
FabArray<FAB>::PostRcvs (const MapOfCopyComTagContainers&  m_RcvVols,
                         const MapOfCopyComTagContainers&  m_RcvTags,
                         Vector<char*>&                     recv_data,
                         Vector<int>&                       recv_size,
                         Vector<int>&                       recv_from,
                         Vector<MPI_Request>&               recv_reqs,
                         int                               icomp,
                         int                               ncomp,
                         int                               SeqNum,
                         int                               preSeqNum,
                         MPI_Comm                          comm)
{
    recv_data.clear();
    recv_size.clear();
    recv_from.clear();
    recv_reqs.clear();

    for (const auto& kv : m_RcvVols) // loop over senders
    {
        std::size_t nbytes = 0;
        if (FAB::preAllocatable())
        {
            for (auto const& cct : kv.second)
            {
                nbytes += (*this)[cct.dstIndex].nBytes(cct.dbox,icomp,ncomp);
            }
        }

        BL_ASSERT(nbytes < std::numeric_limits<int>::max());

        recv_data.push_back(nullptr);
        recv_size.push_back(static_cast<int>(nbytes));
        recv_from.push_back(kv.first);
        recv_reqs.push_back(MPI_REQUEST_NULL);
    }

    const int nrecv = recv_from.size();
    
    Vector<Vector<int> > indv_recv_size;
    Vector<MPI_Request> pre_reqs(nrecv,MPI_REQUEST_NULL);

    if (!FAB::preAllocatable())
    {
        BL_ASSERT(preSeqNum >= 0);
        indv_recv_size.resize(nrecv);
        for (int k = 0; k < nrecv; ++k)
        {
            auto n = m_RcvTags.at(recv_from[k]).size();
            indv_recv_size[k].resize(n);
            pre_reqs[k] = ParallelDescriptor::Arecv(indv_recv_size[k].data(),
                                                    n, recv_from[k], preSeqNum, comm).req();
        }
    }

    int recv_counter = 0;
    while (recv_counter < nrecv)
    {
        int i;

        if (FAB::preAllocatable())
        {
            i = recv_counter;
        }
        else
        {
            MPI_Status status;
            MPI_Waitany(nrecv, pre_reqs.data(), &i, &status);

            for (auto x : indv_recv_size[i])
            {
                recv_size[i] += x;
            }
        }

        if (recv_size[i] > 0)
        {
            recv_data[i] = static_cast<char*>(amrex::The_FA_Arena()->alloc(recv_size[i]));
            recv_reqs[i] = ParallelDescriptor::Arecv(recv_data[i], recv_size[i],
                                                     recv_from[i], SeqNum, comm).req();
        }
        
        ++recv_counter;
    }
}
#endif

#ifdef BL_USE_MPI3
template <class FAB>
void
FabArray<FAB>::PostRcvs_MPI_Onesided (const MapOfCopyComTagContainers&  m_RcvVols,
                                      char*&                            the_recv_data,
                                      Vector<char*>&                     recv_data,
                                      Vector<int>&                       recv_size,
                                      Vector<int>&                       recv_from,
                                      Vector<MPI_Request>&               recv_reqs,
                                      Vector<MPI_Aint> &                 recv_disp,
                                      int                               icomp,
                                      int                               ncomp,
                                      int                               SeqNum,
                                      MPI_Win &                         win)
{
    recv_data.clear();
    recv_size.clear();
    recv_from.clear();
    recv_reqs.clear();
    recv_disp.clear();

    std::size_t TotalRcvsVolume = 0;
    for (const auto& kv : m_RcvVols) // loop over senders
    {
        std::size_t nbytes = 0;
        for (const auto& cct : kv.second)
        {
            nbytes += (*this)[cct.dstIndex].nBytes(cct.dbox,icomp,ncomp);
        }

        BL_ASSERT(nbytes < std::numeric_limits<int>::max());

        TotalRcvsVolume += nbytes;

        recv_data.push_back(nullptr);
        recv_size.push_back(static_cast<int>(nbytes));
        recv_from.push_back(kv.first);
        recv_reqs.push_back(MPI_REQUEST_NULL);
        recv_disp.push_back(0);
    }

    if (TotalRcvsVolume == 0)
    {
        the_recv_data = nullptr;
    }
    else
    {
        the_recv_data = static_cast<char*>(amrex::The_FA_Arena()->alloc(TotalRcvsVolume));

        MPI_Win_attach(win, the_recv_data, TotalRcvsVolume);

        int nrecv = recv_from.size();
        
        std::size_t offset = 0;
        for (int i = 0; i < nrecv; ++i)
        {
            if (recv_size[i] > 0)
            {
                recv_data[i] = the_recv_data + offset;
                MPI_Get_address(recv_data[i], &recv_disp[i]);
                recv_reqs[i] = ParallelDescriptor::Asend(&recv_disp[i],1,recv_from[i],SeqNum).req();
                offset += recv_size[i];
            }
        }
    }
}
#endif   // BL_USE_MPI3

template <class FAB>
bool
FabArray<FAB>::defined (int K) const
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
FabArray<FAB>::defined (const MFIter& mfi) const
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
const FAB&
FabArray<FAB>::operator[] (const MFIter& mfi) const
{
    BL_ASSERT(mfi.LocalIndex() < indexArray.size());
    BL_ASSERT(DistributionMap() == mfi.DistributionMap());
    return *m_fabs_v[mfi.LocalIndex()];
}

template <class FAB>
FAB&
FabArray<FAB>::operator[] (const MFIter& mfi)
{
    BL_ASSERT(mfi.LocalIndex() < indexArray.size());
    BL_ASSERT(DistributionMap() == mfi.DistributionMap());
    return *m_fabs_v[mfi.LocalIndex()];
}

template <class FAB>
const FAB&
FabArray<FAB>::operator[] (int K) const
{
    int li = localindex(K);
    BL_ASSERT(li >=0 && li < indexArray.size());
    return *m_fabs_v[li];
}

template <class FAB>
FAB&
FabArray<FAB>::operator[] (int K)
{
    int li = localindex(K);
    BL_ASSERT(li >=0 && li < indexArray.size());
    return *m_fabs_v[li];
}

template <class FAB>
void
FabArray<FAB>::clear ()
{
    if (define_function_called)
    {
        clearThisBD();  // addThisBD is called in define
    }

    for(Iterator it = m_fabs_v.begin(); it != m_fabs_v.end(); ++it) {
        delete *it;
    }

    m_fabs_v.clear();
    boxarray.clear();
    m_factory.reset();
}

template <class FAB>
template <class>
void
FabArray<FAB>::setVal (value_type val,
                       int        nghost)
{
    setVal(val,0,n_comp,nghost);
}

template <class FAB>
template <class>
void
FabArray<FAB>::setVal (value_type   val,
                       const Box& region,
                       int        nghost)
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
    define(rhs.boxArray(), rhs.DistributionMap(), ncomp, rhs.nGrow(),
           MFInfo().SetAlloc(false), *m_factory);

    if (maketype == amrex::make_alias)
    {
        for (const auto& rhsfab : rhs.m_fabs_v) {
            m_fabs_v.push_back(new FAB(*rhsfab, amrex::make_alias, scomp, ncomp));
        }
    }
    else
    {
        amrex::Abort("FabArray: unknown MakeType");
    }
}

template <class FAB>
template <class>
FabArray<FAB>::FabArray (const BoxArray&            bxs,
                         const DistributionMapping& dm,
                         int                        nvar,
                         int                        ngrow,
                         const Vector<value_type*>&  pval)
    : m_factory(new DefaultFabFactory<FAB>()),
      shmem()
{
    m_FA_stats.recordBuild();
    define(bxs,dm,nvar,ngrow,MFInfo().SetAlloc(false), *m_factory);

    for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
        setFab(mfi, new FAB(mfi.fabbox(),nvar,pval[mfi.LocalIndex()]));
    }
}


template <class FAB>
FabArray<FAB>::FabArray (FabArray<FAB>&& rhs) noexcept
    : FabArrayBase (std::move(rhs))
    , m_factory    (std::move(rhs.m_factory))
    , m_fabs_v     (std::move(rhs.m_fabs_v))
    , shmem        (std::move(rhs.shmem))
    // no need to worry about the data used in non-blocking FillBoundary.
{
    m_FA_stats.recordBuild();
    rhs.m_fabs_v.clear();
    rhs.boxarray.clear();
    rhs.m_factory.reset();
}

template <class FAB>
FabArray<FAB>::~FabArray ()
{
    m_FA_stats.recordDelete();

    if (define_function_called)
    {
	typename std::map<int, std::map<int, FabArray<FAB> *> >::iterator afapIter = 
	    FabArray<FAB>::allocatedFAPointers.find(distributionMap.size());
	if(afapIter == FabArray<FAB>::allocatedFAPointers.end()) {
#ifdef DEBUG_AFAP
	    if(ParallelDescriptor::IOProcessor(this->color()) && distributionMap.size() > 0) {
		std::cout << "**** In FabArray::clear():: map not found:  size = "
			  << distributionMap.size() << std::endl;
	    }
#endif
	} else {
	    std::map<int, FabArray<FAB> *> &faPtrCachedMap = afapIter->second;
	    faPtrCachedMap.erase(aFAPId);
#ifdef DEBUG_AFAP
	    if(ParallelDescriptor::IOProcessor() && faPtrCachedMap.size() == 0) {
		std::cout << "**** In FabArray::~FabArray()::  size distmapngrids = "
			  << faPtrCachedMap.size() << "  " << afapIter-> first << std::endl;
	    }
#endif
	    if(faPtrCachedMap.size() == 0) {  //---- delete unused maps
		FabArray<FAB>::allocatedFAPointers.erase(afapIter);
	    }
	}

	noallocFAPIds.erase(aFAPId);
    }

    clear();
}

template <class FAB>
bool
FabArray<FAB>::ok () const
{
    long isok = true;

    for (MFIter fai(*this); fai.isValid() && isok; ++fai)
    {
        if (defined(fai))
        {
            if (get(fai).box() != fabbox(fai.index()))
            {
                isok = false;
            }
        }
        else
        {
            isok = false;
        }
    }

    ParallelDescriptor::ReduceLongAnd(isok, this->color());

    return isok != 0;
}

template <class FAB>
void
FabArray<FAB>::define (const BoxArray&            bxs,
		       const DistributionMapping& dm,
		       int                        nvar,
		       int                        ngrow,
		       const MFInfo&              info,
                       const FabFactory<FAB>&     factory)
{
    define_function_called = true;

    BL_ASSERT(ngrow >= 0);
    BL_ASSERT(boxarray.size() == 0);
    FabArrayBase::define(bxs, dm, nvar, ngrow);

    addThisBD();

    if(info.alloc) {
        AllocFabs(factory);
#ifdef AMREX_USE_DEVICE
	Device::synchronize();
#endif
    }

    typename std::map<int, std::map<int, FabArray<FAB> *> >::iterator afapIter = 
                 FabArray<FAB>::allocatedFAPointers.find(distributionMap.size());
    if(afapIter == FabArray<FAB>::allocatedFAPointers.end()) {
      int dmapSize(distributionMap.size());
      std::map<int, FabArray<FAB> *> tempMap;
      tempMap.insert(std::make_pair(aFAPId, this));
      FabArray<FAB>::allocatedFAPointers.insert(std::make_pair(dmapSize, tempMap));

    } else {
      std::map<int, FabArray<FAB> *> &faPtrCachedMap = afapIter->second;
      faPtrCachedMap.insert(std::make_pair(aFAPId, this));
    }

    if(info.alloc) {
      noallocFAPIds.erase(aFAPId);
    } else {
      noallocFAPIds.insert(aFAPId);
    }

#ifdef BL_USE_TEAM
    if (info.alloc) ParallelDescriptor::MyTeam().MemoryBarrier();
#endif
}

template <class FAB>
void
FabArray<FAB>::AllocFabs (const FabFactory<FAB>& factory)
{
    const int n = indexArray.size();
    const int nworkers = ParallelDescriptor::TeamSize();
    shmem.alloc = (nworkers > 1);

    bool alloc = !shmem.alloc;

    FabInfo fab_info;
    fab_info.SetAlloc(alloc).SetShared(shmem.alloc);

    m_factory.reset(factory.clone());    

    m_fabs_v.reserve(n);

    for (int i = 0; i < n; ++i)
    {
	int K = indexArray[i];
        const Box& tmpbox = fabbox(K);
        m_fabs_v.push_back(m_factory->create(tmpbox, n_comp, fab_info, K));
    }
    
#ifdef BL_USE_TEAM
    if (shmem.alloc)
    {
	const int teamlead = ParallelDescriptor::MyTeamLead();

	shmem.n_values = 0;
	shmem.n_points = 0;
	Vector<long> offset(n,0);
	Vector<long> nextoffset(nworkers,-1);
	for (int i = 0; i < n; ++i) {
	    int K = indexArray[i];
	    int owner = distributionMap[K] - teamlead;
	    long s = m_fabs_v[i]->size();
	    if (ownership[i]) {
		shmem.n_values += s;
		shmem.n_points += m_fabs_v[i]->nPts();
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

#if defined (BL_USE_UPCXX)

	shmem.p = BLPgas::alloc(bytes);

	upcxx::global_ptr<void> psrc = upcxx::global_ptr<void>(shmem.p);
	std::vector<upcxx::global_ptr<void> > pdst(nworkers);
	const auto& team = ParallelDescriptor::MyTeam().get();
	team.allgather(&psrc, &pdst[0], sizeof(upcxx::global_ptr<void>));

	mfp = static_cast<value_type*>(shmem.p);

	for (int w = 0; w < nworkers; ++w) {
	    void* p = (void*) pdst[w];
	    dps.push_back((value_type*)p);
	}

#elif defined (BL_USE_MPI3)

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

	amrex::Abort("BaseFab::define: to allocate shared memory, either USE_UPCXX or USE_MPI3 must be true");

#endif

	for (int i = 0; i < n; ++i) {
	    int K = indexArray[i];
	    int owner = distributionMap[K] - teamlead;
	    value_type *p = dps[owner] + offset[i];
	    m_fabs_v[i]->setPtr(p, m_fabs_v[i]->size());
	}

	for (long i = 0; i < shmem.n_values; i++, mfp++) {
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
	m_fabs_v.resize(indexArray.size());
    }

    m_fabs_v[localindex(boxno)] = elem;
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
        m_fabs_v.resize(indexArray.size());
    }

    m_fabs_v[mfi.LocalIndex()] = elem;
}

template <class FAB>
template <class>
void
FabArray<FAB>::setBndry (value_type val)
{
    setBndry(val, 0, n_comp);
}

template <class FAB>
template <class>
void
FabArray<FAB>::setBndry (value_type val,
                         int        strt_comp,
                         int        ncomp)
{
    if (n_grow > 0)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter fai(*this); fai.isValid(); ++fai)
        {
            get(fai).setComplement(val, fai.validbox(), strt_comp, ncomp);
        }
    }
}

template <class FAB>
template <class>
void
FabArray<FAB>::setDomainBndry (value_type val, const Geometry& geom)
{
    setDomainBndry(val, 0, n_comp, geom);
}

template <class FAB>
template <class>
void
FabArray<FAB>::setDomainBndry (value_type val,
                               int        strt_comp,
                               int        ncomp,
                               const Geometry& geom)
{
    if (n_grow > 0)
    {
        const Box& domain_box = amrex::convert(geom.Domain(), boxArray().ixType());
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter fai(*this); fai.isValid(); ++fai)
        {
            get(fai).setComplement(val, domain_box, strt_comp, ncomp);
        }
    }
}

template <class FAB>
void
FabArray<FAB>::ParallelCopy (const FabArray<FAB>& src,
                             int                  scomp,
                             int                  dcomp,
                             int                  ncomp,
                             int                  snghost,
                             int                  dnghost,
                             const Periodicity&   period,
                             CpOp                 op)
{
    BL_PROFILE("FabArray::ParallelCopy()");

    if (size() == 0 || src.size() == 0) return;

    BL_ASSERT(op == FabArrayBase::COPY || op == FabArrayBase::ADD);
    BL_ASSERT(boxArray().ixType() == src.boxArray().ixType());

    BL_ASSERT(src.nGrow() >= snghost);
    BL_ASSERT(    nGrow() >= dnghost);

    if ((src.boxArray().ixType().cellCentered() || op == FabArrayBase::COPY) &&
        (boxarray == src.boxarray && distributionMap == src.distributionMap)
	&& snghost ==0 && dnghost == 0 && !period.isAnyPeriodic())
    {
        //
        // Short-circuit full intersection code if we're doing copy()s or if
        // we're doing plus()s on cell-centered data.  Don't do plus()s on
        // non-cell-centered data this simplistic way.
        //
#ifdef _OPENMP
#pragma omp parallel if (FAB::isCopyOMPSafe())
#endif
        for (MFIter fai(*this,true); fai.isValid(); ++fai)
        {
            const Box& bx = fai.tilebox();

	    if (this != &src) {
		// avoid self copy or plus
		if (op == FabArrayBase::COPY) {
		    get(fai).copy(src[fai],bx,scomp,bx,dcomp,ncomp);
		} else {
		    get(fai).plus(src[fai],bx,bx,scomp,dcomp,ncomp);
		}
	    }
        }

        return;
    }

    const CPC& thecpc = getCPC(dnghost, src, snghost, period);

    if (ParallelDescriptor::NProcs() == 1)
    {
        //
        // There can only be local work to do.
        //
	int N_loc = (*thecpc.m_LocTags).size();
#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe() && thecpc.m_threadsafe_loc)
#endif
	for (int i=0; i<N_loc; ++i)
        {
            const CopyComTag& tag = (*thecpc.m_LocTags)[i];

	    if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
		// avoid self copy or plus
		if (op == FabArrayBase::COPY) {
		    get(tag.dstIndex).copy(src[tag.srcIndex],tag.sbox,scomp,tag.dbox,dcomp,ncomp);
		} else {
		    get(tag.dstIndex).plus(src[tag.srcIndex],tag.sbox,tag.dbox,scomp,dcomp,ncomp);
		}
	    }
        }

#ifdef AMREX_USE_DEVICE
        Device::synchronize();
#endif

        return;
    }

#ifdef BL_USE_MPI

#if defined(BL_USE_UPCXX)
    ParallelDescriptor::Mode.set_upcxx_mode();
#endif

#if defined(BL_USE_MPI3)
    BL_ASSERT(FAB::preAllocatable());
#else
    BL_ASSERT(!ParallelDescriptor::MPIOneSided());
#endif

    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    int SeqNum, preSeqNum;
    {
	ParallelDescriptor::Color src_color = src.color();
	ParallelDescriptor::Color dst_color = this->color();
	if (src_color == ParallelDescriptor::DefaultColor() ||
	    dst_color == ParallelDescriptor::DefaultColor() ||
	    src_color != dst_color) {
	    // If any of the two FabArrays is in the global communicator, 
	    // or if the two have different color,
	    // all processes are here.
            if (!FAB::preAllocatable()) {
                preSeqNum = ParallelDescriptor::SeqNum();
            }
	    SeqNum  = ParallelDescriptor::SeqNum();
	} else { // The two have the same non-default color.
	    if (ParallelDescriptor::SubCommColor() == src_color) {
                if (!FAB::preAllocatable()) {
                    preSeqNum = ParallelDescriptor::SubSeqNum();
                }
		SeqNum  = ParallelDescriptor::SubSeqNum();
	    }
	    // else I don't have any data and my SubSeqNum() should not be called.
	}
    }	

    const int N_snds = thecpc.m_SndTags->size();
    const int N_rcvs = thecpc.m_RcvTags->size();
    const int N_locs = thecpc.m_LocTags->size();

    if (N_locs == 0 && N_rcvs == 0 && N_snds == 0)
        //
        // No work to do.
        //
        return;

#ifdef BL_USE_MPI3
    MPI_Group tgroup, rgroup, sgroup;
    if (ParallelDescriptor::MPIOneSided()) {
	MPI_Comm_group(ParallelDescriptor::Communicator(this->color()), &tgroup);
    }
#endif

    //
    // Send/Recv at most MaxComp components at a time to cut down memory usage.
    //
    int NCompLeft = ncomp;

    for (int ipass = 0, SC = scomp, DC = dcomp; ipass < ncomp; )
    {
        const int NC = std::min(NCompLeft,FabArrayBase::MaxComp);

        //
        // Before we post recv, let's preprocess sends in case FAB is not preAllocatable
        //

	Vector<char*>                       send_data;
	Vector<int>                         send_size;
	Vector<int>                         send_rank;
	Vector<MPI_Request>                 send_reqs;
	Vector<const CopyComTagsContainer*> send_cctc;
        Vector<Vector<int> >                 indv_send_size;
        Vector<MPI_Request>                 pre_reqs;

#if defined BL_USE_UPCXX || defined BL_USE_MPI3 
        int actual_n_snds = 0;
#endif
	if (N_snds > 0)
	{
	    send_data.reserve(N_snds);
	    send_size.reserve(N_snds);
	    send_rank.reserve(N_snds);
            send_reqs.reserve(N_snds);
	    send_cctc.reserve(N_snds);
            indv_send_size.reserve(N_snds);

            for (auto const& kv : *thecpc.m_SndVols)
	    {
                Vector<int> iss;                
                auto const& cctc = thecpc.m_SndTags->at(kv.first);

                std::size_t nbytes = 0;
                if (FAB::preAllocatable())
                {
                    for (auto const& cct : kv.second)
                    {
                        nbytes += src[cct.srcIndex].nBytes(cct.sbox,SC,NC);
                    }
                }
                else
                {
                    for (auto const& tag : cctc)
                    {
                        std::size_t b = src[tag.srcIndex].nBytes(tag.sbox,SC,NC);
                        nbytes += b;
                        iss.push_back(static_cast<int>(b));
                    }
                }
		
		BL_ASSERT(nbytes < std::numeric_limits<int>::max());

                char* data = nullptr;
                if (nbytes > 0)
                {
                    data = static_cast<char*>
#ifdef BL_USE_UPCXX
                        (BLPgas::alloc(nbytes));
#else
                        (amrex::The_FA_Arena()->alloc(nbytes));
#endif
                }
                    
		send_data.push_back(data);
                send_size.push_back(static_cast<int>(nbytes));
                send_rank.push_back(kv.first);
                send_reqs.push_back(MPI_REQUEST_NULL);
                send_cctc.push_back(&cctc);
                indv_send_size.push_back(std::move(iss));
	    }

#if defined BL_USE_UPCXX || defined BL_USE_MPI3
            actual_n_snds = N_snds - std::count(send_size.begin(), send_size.end(), 0);
#endif
        }

        if (!FAB::preAllocatable())
        {
            pre_reqs.resize(N_snds,MPI_REQUEST_NULL);
            for (int j=0; j<N_snds; ++j)
            {
                pre_reqs[j] = ParallelDescriptor::Asend(indv_send_size[j].data(),
                                                        indv_send_size[j].size(),
                                                        send_rank[j],preSeqNum).req();
            }
        }

        Vector<int>         recv_from;
        Vector<char*>       recv_data;
        Vector<int>         recv_size;
        Vector<MPI_Request> recv_reqs;
#ifdef BL_USE_MPI3
	Vector<MPI_Aint>    recv_disp;
#endif
        //
        // Post rcvs. Allocate one chunk of space to hold'm all.
        //
        char* the_recv_data = nullptr;

        int actual_n_rcvs = 0;
	if (N_rcvs > 0) {
#ifdef BL_USE_UPCXX
	    PostRcvs_PGAS(*thecpc.m_RcvVols, the_recv_data, recv_data,
                          recv_size, recv_from, 
                          SC, NC, SeqNum, &BLPgas::cp_recv_event);
#else
	    if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
                PostRcvs_MPI_Onesided(*thecpc.m_RcvVols, the_recv_data, recv_data, 
                                      recv_size, recv_from, recv_reqs, recv_disp,
                                      SC, NC, SeqNum, ParallelDescriptor::cp_win);
		MPI_Group_incl(tgroup, recv_from.size(), recv_from.dataPtr(), &rgroup);
		MPI_Win_post(rgroup, 0, ParallelDescriptor::cp_win);
#endif
	    } else {
                PostRcvs(*thecpc.m_RcvVols, *thecpc.m_RcvTags,
                         recv_data, recv_size, recv_from, recv_reqs, SC, NC, SeqNum, preSeqNum);
	    }
#endif
            actual_n_rcvs = N_rcvs - std::count(recv_size.begin(), recv_size.end(), 0);
	}

	//
	// Post send's
	// 
	if (N_snds > 0)
	{
#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe())
#endif
	    for (int j=0; j<N_snds; ++j)
	    {
                char* dptr = send_data[j];
                if (dptr != nullptr)
                {
                    auto const& cctc = *send_cctc[j];
                    for (auto const& tag : cctc)
                    {
                        const Box& bx = tag.sbox;
                        auto n = src[tag.srcIndex].copyToMem(bx,SC,NC,dptr);
                        dptr += n;
                    }
                    BL_ASSERT(dptr == send_data[j] + send_size[j]);
		}
	    }

#ifdef AMREX_USE_DEVICE
            Device::synchronize();
#endif

#ifdef BL_USE_UPCXX
	    
	    BLPgas::cp_send_counter = 0;

	    for (int i=0; i<N_snds; ++i) {
                if (send_size[i] > 0) {
                    BLPgas::Send(upcxx::global_ptr<void>((void *)send_data[i], upcxx::myrank()),
                                 send_rank[i], send_size[i], SeqNum,
                                 &BLPgas::cp_send_event, &BLPgas::cp_send_counter);
                }
	    }
	    
	    // Need to make sure at least half of the sends have been started
	    while (BLPgas::cp_send_counter < actual_n_snds) {
		upcxx::advance();
            }

#else
	    
	    if (ParallelDescriptor::MPIOneSided())
	    {
#if defined(BL_USE_MPI3)
		Vector<MPI_Aint> send_disp(N_snds,0);
                send_reqs.clear();

		for (int i=0; i<N_snds; ++i) {
                    if (send_size[i] > 0) {
                        send_reqs[i] = ParallelDescriptor::Arecv
                            (&send_disp[i],1,send_rank[i],SeqNum).req();
                    }
		}
		
		MPI_Group_incl(tgroup, N_snds, send_rank.data(), &sgroup);
		MPI_Win_start(sgroup,0,ParallelDescriptor::cp_win);
		
		int send_counter = 0;	
		while (send_counter < actual_n_snds) {
		    MPI_Status status;
		    int index;
		    
		    MPI_Waitany(N_snds, send_reqs.dataPtr(), &index, &status);
		    
		    BL_ASSERT(status.MPI_TAG == SeqNum);
		    BL_ASSERT(status.MPI_SOURCE == send_rank[index]);
		    
		    MPI_Put(send_data[index], send_size[index], MPI_CHAR, send_rank[index],
			    send_disp[index], send_size[index], MPI_CHAR, ParallelDescriptor::cp_win);
		    
		    ++send_counter;
		}
#endif
	    }
            else
            {
                int send_counter = 0;
                while (send_counter < N_snds)
                {
                    int j;

                    if (FAB::preAllocatable())
                    {
                        j = send_counter;
                    }
                    else
                    {
                        MPI_Status status;
                        MPI_Waitany(N_snds, pre_reqs.data(), &j, &status);
                    }
                        
                    if (send_size[j] > 0)
                    {
                        if (FabArrayBase::do_async_sends)
                        {
                            send_reqs[j] = ParallelDescriptor::Asend
                                (send_data[j],send_size[j],send_rank[j],SeqNum).req();
                        }
                        else
                        {
                            ParallelDescriptor::Send(send_data[j],send_size[j],send_rank[j],SeqNum);
                            amrex::The_FA_Arena()->free(send_data[j]);
                        }
                    }

                    ++send_counter;
		}
	    }
#endif
	}

#ifdef BL_USE_MPI3
	if (ParallelDescriptor::MPIOneSided()) {
	    if (N_rcvs > 0) MPI_Group_free(&rgroup);
	    if (N_snds > 0) MPI_Group_free(&sgroup);
	}
#endif

        //
        // Do the local work.  Hope for a bit of communication/computation overlap.
        //
	if (ParallelDescriptor::TeamSize() > 1 && thecpc.m_threadsafe_loc)
	{
#ifdef BL_USE_TEAM
#ifdef _OPENMP
#pragma omp parallel if (FAB::isCopyOMPSafe())
#endif
	    ParallelDescriptor::team_for(0, N_locs, [&] (int j) 
            {
		const CopyComTag& tag = (*thecpc.m_LocTags)[j];
		
		if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
		    // avoid self copy or plus
		    if (op == FabArrayBase::COPY) {
			get(tag.dstIndex).copy(src[tag.srcIndex],tag.sbox,SC,tag.dbox,DC,NC);
		    } else {
			get(tag.dstIndex).plus(src[tag.srcIndex],tag.sbox,tag.dbox,SC,DC,NC);
		    }
		}
	    });
#endif	    
	}
	else 
	{
#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe() && thecpc.m_threadsafe_loc)
#endif
	    for (int j=0; j<N_locs; ++j)
	    {
		const CopyComTag& tag = (*thecpc.m_LocTags)[j];

		if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
		    // avoid self copy or plus
		    if (op == FabArrayBase::COPY) {
			get(tag.dstIndex).copy(src[tag.srcIndex],tag.sbox,SC,tag.dbox,DC,NC);
		    } else {
			get(tag.dstIndex).plus(src[tag.srcIndex],tag.sbox,tag.dbox,SC,DC,NC);
		    }
		}
	    }
	}

#ifdef AMREX_USE_DEVICE
        Device::synchronize();
#endif

	//
	//  wait and unpack
	//
#ifdef BL_USE_UPCXX
        if (actual_n_rcvs > 0) BLPgas::cp_recv_event.wait();
#else
	if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
	    if (N_snds > 0) MPI_Win_complete(ParallelDescriptor::cp_win);
	    if (N_rcvs > 0) MPI_Win_wait    (ParallelDescriptor::cp_win);
#endif
	} else {
	    if (actual_n_rcvs > 0) {
		Vector<MPI_Status> stats(N_rcvs);
		BL_MPI_REQUIRE( MPI_Waitall(N_rcvs, recv_reqs.dataPtr(), stats.dataPtr()) );
		if (!CheckRcvStats(stats, recv_size, MPI_CHAR, SeqNum))
                {
                    amrex::Abort("ParallelCopy failed with wrong message size");
                }
	    }
	}
#endif	

	if (N_rcvs > 0)
	{
	    Vector<const CopyComTagsContainer*> recv_cctc(N_rcvs,nullptr);

	    for (int k = 0; k < N_rcvs; ++k)
	    {
                if (recv_size[k] > 0)
                {
                    auto const& cctc = thecpc.m_RcvTags->at(recv_from[k]);
                    recv_cctc[k] = &cctc;
                }
	    }
	    
#ifdef _OPENMP
#pragma omp parallel if (FAB::isCopyOMPSafe() && thecpc.m_threadsafe_rcv)
#endif
	    {
		FAB fab;

#ifdef _OPENMP
#pragma omp for
#endif
                for (int k = 0; k < N_rcvs; ++k)
		{
		    const char* dptr = recv_data[k];
                    if (dptr != nullptr)
                    {
                        auto const& cctc = *recv_cctc[k];
                        for (auto const& tag : cctc)
                        {
                            const Box& bx  = tag.dbox;
                            std::size_t n;
                            if (op == FabArrayBase::COPY)
                            {
                                n = get(tag.dstIndex).copyFromMem(bx,DC,NC,dptr);
                            }
                            else
                            {
                                fab.resize(bx,NC);
                                n = fab.copyFromMem(bx,0,NC,dptr);
                                get(tag.dstIndex).plus(fab,bx,bx,0,DC,NC);
                            }
                            dptr += n;
                        }
                        BL_ASSERT(dptr == recv_data[k] + recv_size[k]);
		    }
		}
	    }

#ifdef AMREX_USE_DEVICE
            Device::synchronize();
#endif

            if (the_recv_data)
            {
#ifdef BL_USE_UPCXX
                BLPgas::free(the_recv_data);
#else
                if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
                    MPI_Win_detach(ParallelDescriptor::cp_win, the_recv_data);
#endif
                }
                amrex::The_FA_Arena()->free(the_recv_data);
#endif
            }
            else
            {
                for (auto p : recv_data) {
                    amrex::The_FA_Arena()->free(p);
                }
            }
	}
	
        if (N_snds > 0) {
#ifdef  BL_USE_UPCXX
	    FabArrayBase::WaitForAsyncSends_PGAS(N_snds,send_data,
					         &BLPgas::cp_send_event,
					         &BLPgas::cp_send_counter);
#else
	    if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
		for (int i = 0; i < N_snds; ++i) {
		    if (send_data[i]) amrex::The_FA_Arena()->free(send_data[i]);
                }
#endif
	    } else {
		if (FabArrayBase::do_async_sends && ! thecpc.m_SndTags->empty()) {
		    Vector<MPI_Status> stats;
		    FabArrayBase::WaitForAsyncSends(N_snds,send_reqs,send_data,stats);
		}
	    }
#endif
        }

        ipass     += NC;
        SC        += NC;
        DC        += NC;
        NCompLeft -= NC;
    }

#ifdef BL_USE_MPI3
    if (ParallelDescriptor::MPIOneSided()) {
	MPI_Group_free(&tgroup);
    }
#endif

#ifdef BL_USE_TEAM
    ParallelDescriptor::MyTeam().MemoryBarrier();
#endif

    return;

#endif /*BL_USE_MPI*/
}

template <class FAB>
void
FabArray<FAB>::ParallelCopy (const FabArray<FAB>& src,
                             int                  scomp,
                             int                  dcomp,
                             int                  ncomp,
                             const Periodicity&   period,
                             CpOp                 op)
{
    ParallelCopy(src,scomp,dcomp,ncomp,0,0,period,op);
}

template <class FAB>
void
FabArray<FAB>::ParallelCopy (const FabArray<FAB>& src, const Periodicity& period, CpOp op)
{
    copy(src,0,0,nComp(),0,0,period,op);
}

//
// Copies to FABs, note that destination is first arg.
//

template <class FAB>
void
FabArray<FAB>::copyTo (FAB& dest,
		       int  nghsot) const
{
    copyTo(dest, dest.box(), 0, 0, dest.nComp(), nghsot);
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
void
FabArray<FAB>::copyTo (FAB&       dest,
		       const Box& subbox,
		       int        scomp,
		       int        dcomp,
		       int        ncomp,
		       int        nghost) const
{
    BL_PROFILE("FabArray::copy(fab)");

    BL_ASSERT(dcomp + ncomp <= dest.nComp());
    BL_ASSERT(nghost <= nComp());

    if (ParallelDescriptor::NProcs() == 1)
    {
        for (int j = 0, N = size(); j < N; ++j)
        {
	    const Box& bx = amrex::grow(boxarray[j],nghost);
	    const Box& destbox = bx & subbox;
	    if (destbox.ok())
            {
                dest.copy(get(j),destbox,scomp,destbox,dcomp,ncomp);
            }
        }

        return;
    }

    //
    //  Note that subbox must be identical on each process!!
    //
#ifdef DEBUG
    {
	BoxCommHelper bch(subbox);	
	BL_ASSERT(this->color() == ParallelDescriptor::DefaultColor());
	ParallelDescriptor::Bcast(bch.data(), bch.size(), 0);
	const Box& bx0 = bch.make_box();
	BL_ASSERT(subbox == bx0);
    }
#endif

    FAB ovlp;

    std::vector< std::pair<int,Box> > isects;
    boxarray.intersections(subbox, isects, false, nghost);

    for (int j = 0, M = isects.size(); j < M; ++j)
    {
	const int  k  = isects[j].first;
	const Box& bx = isects[j].second;

	ovlp.resize(bx,ncomp);

	if (ParallelDescriptor::MyProc() == distributionMap[k])
	{
	    ovlp.copy(get(k),bx,scomp,bx,0,ncomp);
	}

	const int N = bx.numPts()*ncomp;

	BL_ASSERT(this->color() == ParallelDescriptor::DefaultColor());
	ParallelDescriptor::Bcast(ovlp.dataPtr(),N,distributionMap[k]);

	dest.copy(ovlp,bx,0,bx,dcomp,ncomp);
    }
}




template <class FAB>
void
FabArray<FAB>::copyInter (FabArray<FAB> *src, FabArray<FAB> *dest,
                          int                  scomp,
                          int                  dcomp,
                          int                  ncomp,
                          int                  snghost,
                          int                  dnghost,
                          const MPI_Comm      &commSrc,
                          const MPI_Comm      &commDest,
                          const MPI_Comm      &commInter,
                          const MPI_Comm      &commBoth,
                	  bool                 isSrc,
                          CpOp                 op)
{
#ifdef BL_USE_MPI
    BL_PROFILE("FabArray::copyInter()");

    bool isDest( ! isSrc);
    if(isSrc) {
      BL_ASSERT(src  != 0);
      BL_ASSERT(dest == 0);
    }
    if(isDest) {
      BL_ASSERT(src  == 0);
      BL_ASSERT(dest != 0);
    }

    if ((isDest && dest->size() == 0) || (isSrc && src->size() == 0)) return;

    BL_ASSERT(op == FabArrayBase::COPY || op == FabArrayBase::ADD);

    if(isSrc)  { BL_ASSERT(src->nGrow() >= snghost);  }
    if(isDest) { BL_ASSERT(dest->nGrow() >= dnghost); }

    MPI_Group groupSrc(MPI_GROUP_NULL), groupDest(MPI_GROUP_NULL), groupBoth(MPI_GROUP_NULL);
    MPI_Group groupAll(MPI_GROUP_NULL);
    int myProcAll(ParallelDescriptor::MyProcAll()), myProcBoth(MPI_UNDEFINED);
    int iopN(ParallelDescriptor::IOProcessorNumber()), iopNInBoth(MPI_UNDEFINED);
    int bcastRootSrc(MPI_UNDEFINED), bcastRootDest(MPI_UNDEFINED);
    size_t sizeOne(1);
    int tagOne(1), tagTwo(2);

    BL_MPI_REQUIRE( MPI_Comm_group(ParallelDescriptor::CommunicatorAll(), &groupAll) );
    BL_MPI_REQUIRE( MPI_Comm_group(commBoth, &groupBoth) );

    BL_MPI_REQUIRE( MPI_Group_translate_ranks(groupAll, sizeOne, &myProcAll, groupBoth, &myProcBoth) );

    if(isSrc) {
      BL_MPI_REQUIRE( MPI_Comm_group(commSrc, &groupSrc) );
      BL_MPI_REQUIRE( MPI_Group_translate_ranks(groupSrc, sizeOne, &iopN, groupBoth, &iopNInBoth) );
      if(ParallelDescriptor::IOProcessor()) {
        ParallelDescriptor::Asend(&iopNInBoth, sizeOne, iopN, tagOne, commBoth);
      }
    }
    if(isDest) {
      BL_MPI_REQUIRE( MPI_Comm_group(commDest, &groupDest) );
      BL_MPI_REQUIRE( MPI_Group_translate_ranks(groupDest, sizeOne, &iopN, groupBoth, &iopNInBoth) );
      if(ParallelDescriptor::IOProcessor()) {
        ParallelDescriptor::Asend(&iopNInBoth, sizeOne, iopN, tagTwo, commBoth);
      }
    }
    if(myProcBoth == iopN) {
      ParallelDescriptor::Recv(&bcastRootSrc,  sizeOne, MPI_ANY_SOURCE, tagOne, commBoth);
      ParallelDescriptor::Recv(&bcastRootDest, sizeOne, MPI_ANY_SOURCE, tagTwo, commBoth);
    }
    ParallelDescriptor::Bcast(&bcastRootSrc,  sizeOne, iopN, commBoth);
    ParallelDescriptor::Bcast(&bcastRootDest, sizeOne, iopN, commBoth);


    // ---- give all procs in commBoth source and dest boxarrays and distribution maps
    BoxArray srcBoxArray, destBoxArray;
    Vector<int> srcPMap, destPMap, srcPMapAll, destPMapAll;
    DistributionMapping srcDM, destDM, srcDMAll, destDMAll;
    Vector<int> srcIndexArray, destIndexArray;
    if(isSrc) {
      srcIndexArray = src->IndexArray();
    }
    if(isDest) {
      destIndexArray = dest->IndexArray();
    }


    if(isSrc && myProcBoth == bcastRootSrc) {
      srcBoxArray = src->boxArray();
      srcDM = src->DistributionMap();
      srcPMap = srcDM.ProcessorMap();
      srcPMapAll = DistributionMapping::TranslateProcMap(srcPMap, groupBoth, groupSrc);
    }
    if(isDest && myProcBoth == bcastRootDest) {
      destBoxArray = dest->boxArray();
      destDM = dest->DistributionMap();
      destPMap = destDM.ProcessorMap();
      destPMapAll = DistributionMapping::TranslateProcMap(destPMap, groupBoth, groupDest);
    }

    amrex::BroadcastBoxArray(srcBoxArray,  myProcBoth, bcastRootSrc,  commBoth);

    amrex::BroadcastBoxArray(destBoxArray, myProcBoth, bcastRootDest, commBoth);
    BL_ASSERT(destBoxArray[0].ixType() == srcBoxArray[0].ixType());


    amrex::BroadcastArray(srcPMapAll,  myProcBoth, bcastRootSrc,  commBoth);
    amrex::BroadcastArray(destPMapAll, myProcBoth, bcastRootDest, commBoth);

    srcDMAll.define(srcPMapAll);
    destDMAll.define(destPMapAll);

    CPC thecpc(destBoxArray, destDMAll, destIndexArray, dnghost,
	       srcBoxArray ,  srcDMAll,  srcIndexArray, dnghost, 
	       Periodicity::NonPeriodic(), myProcBoth);

    if (ParallelDescriptor::NProcsAll() == 1) {
        return;
    }

    // ---- syncronize sequence numbers across groups
    int currentSeqNumber(-4);
    if(isSrc && myProcBoth == bcastRootSrc) {
      currentSeqNumber = ParallelDescriptor::SeqNum(1);
    }
    ParallelDescriptor::Bcast(&currentSeqNumber, sizeOne, bcastRootSrc, commBoth);
    ParallelDescriptor::SeqNum(2, currentSeqNumber);

    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    int preSeqNum;
    if (!FAB::preAllocatable()) {
        preSeqNum = ParallelDescriptor::SeqNum();
    }
    const int SeqNum = ParallelDescriptor::SeqNum();

    if (thecpc.m_LocTags->empty() && thecpc.m_RcvTags->empty() && thecpc.m_SndTags->empty()) {
        //
        // No work to do.
        //
        return;
    }
    //
    // Send/Recv at most MaxComp components at a time to cut down memory usage.
    //
    int NCompLeft = ncomp;

    for (int ipass = 0, SC = scomp, DC = dcomp; ipass < ncomp; )
    {
        const int NC = std::min(NCompLeft,FabArrayBase::MaxComp);

        //
        // Before we post recv, let's preprocess sends in case FAB is not preAllocatable
        //
	const int N_snds = thecpc.m_SndTags->size();

	Vector<char*>                       send_data;
	Vector<int>                         send_size;
	Vector<int>                         send_rank;
	Vector<MPI_Request>                 send_reqs;
	Vector<const CopyComTagsContainer*> send_cctc;
        Vector<Vector<int> >                 indv_send_size;
        Vector<MPI_Request>                 pre_reqs;

#if defined BL_USE_UPCXX || defined BL_USE_MPI3
        int actual_n_snds = 0;
#endif
        if (N_snds > 0)
        {
            send_data.reserve(N_snds);
            send_size.reserve(N_snds);
            send_rank.reserve(N_snds);
            send_reqs.reserve(N_snds);
            send_cctc.reserve(N_snds);
            indv_send_size.reserve(N_snds);
            
            for (auto const& kv : *thecpc.m_SndVols)
            {
                Vector<int> iss;                
                auto const& cctc = thecpc.m_SndTags->at(kv.first);

                std::size_t nbytes = 0;
                if (FAB::preAllocatable())
                {
                    for (auto const& cct : kv.second)
                    {
                        nbytes += (*src)[cct.srcIndex].nBytes(cct.sbox,SC,NC);
                    }
                }
                else
                {
                    for (auto const& tag : cctc)
                    {
                        std::size_t b = (*src)[tag.srcIndex].nBytes(tag.sbox,SC,NC);
                        nbytes += b;
                        iss.push_back(static_cast<int>(b));
                    }
                }
                
                BL_ASSERT(nbytes < std::numeric_limits<int>::max());
                
                char* data = nullptr;
                if (nbytes > 0) {
                    data = static_cast<char*>(amrex::The_FA_Arena()->alloc(nbytes));
                }
                    
                send_data.push_back(data);
                send_size.push_back(static_cast<int>(nbytes));
                send_rank.push_back(kv.first);
                send_reqs.push_back(MPI_REQUEST_NULL);
                send_cctc.push_back(&cctc);
                indv_send_size.push_back(std::move(iss));
            }

#if defined BL_USE_UPCXX || defined BL_USE_MPI3
            actual_n_snds = N_snds - std::count(send_size.begin(), send_size.end(), 0);
#endif
        }

        if (!FAB::preAllocatable())
        {
            pre_reqs.resize(N_snds,MPI_REQUEST_NULL);
            for (int j=0; j<N_snds; ++j)
            {
                pre_reqs[j] = ParallelDescriptor::Asend(indv_send_size[j].data(),
                                                        indv_send_size[j].size(),
                                                        send_rank[j],preSeqNum).req();
            }
        }

        Vector<MPI_Status>  stats;
        Vector<int>         recv_from;
        Vector<char*>       recv_data;
        Vector<int>         recv_size;
        Vector<MPI_Request> recv_reqs;
        //
        // Post rcvs. Allocate one chunk of space to hold'm all.
        //
        char* the_recv_data = nullptr;
        const int N_rcvs = thecpc.m_RcvTags->size();
        int actual_n_rcvs = 0;
        if (N_rcvs > 0) {
            dest->PostRcvs(*thecpc.m_RcvVols, *thecpc.m_RcvTags,
                           recv_data, recv_size, recv_from, recv_reqs, 
                           SC, NC, SeqNum, preSeqNum, commBoth);
            actual_n_rcvs = N_rcvs - std::count(recv_size.begin(), recv_size.end(), 0);
        }

	//
	// Post sends
	// 
        if (N_snds > 0)
        {
#ifdef _OPENMP
#pragma omp parallel for if(FAB::isCopyOMPSafe())
#endif
            for (int j=0; j<N_snds; ++j)
            {
                char* dptr = send_data[j];
                if (dptr != nullptr)
                {
                    auto const& cctc = *send_cctc[j];
                    for (auto const& tag : cctc)
                    {
                        const Box& bx = tag.sbox;
                        auto n = (*src)[tag.srcIndex].copyToMem(bx,SC,NC,dptr);
                        dptr += n;
                    }
                    BL_ASSERT(dptr == send_data[j] + send_size[j]);
                }
            }

#ifdef AMREX_USE_DEVICE
            Device::synchronize();
#endif

            int send_counter = 0;
            while (send_counter < N_snds)
            {
                int j;

                if (FAB::preAllocatable())
                {
                    j = send_counter;
                }
                else
                {
                    MPI_Status status;
                    MPI_Waitany(N_snds, pre_reqs.data(), &j, &status);
                }
                
                if (send_size[j] > 0)
                {
                    if (FabArrayBase::do_async_sends)
                    {
                        send_reqs[j] = ParallelDescriptor::Asend
                            (send_data[j],send_size[j],send_rank[j],SeqNum, commBoth).req();
                    }
                    else
                    {
                        ParallelDescriptor::Send(send_data[j],send_size[j],send_rank[j],SeqNum, commBoth);
                        amrex::The_FA_Arena()->free(send_data[j]);
                    }
                }

                ++send_counter;
            }
        }

        //
        // Do the local work.  Hope for a bit of communication/computation overlap.
        //
	int N_loc = (*thecpc.m_LocTags).size();

#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe() && thecpc.m_threadsafe_loc)
#endif
	for (int j=0; j<N_loc; ++j)
        {
            const CopyComTag& tag = (*thecpc.m_LocTags)[j];

            if (op == FabArrayBase::COPY)
            {
                dest->get(tag.dstIndex).copy((*src)[tag.srcIndex],tag.sbox,SC,tag.dbox,DC,NC);
            }
            else
            {
                dest->get(tag.dstIndex).plus((*src)[tag.srcIndex],tag.sbox,tag.dbox,SC,DC,NC);
            }
        }

	//
	//  wait and unpack
	//

	if (N_rcvs > 0)
	{
	    Vector<const CopyComTagsContainer*> recv_cctc(N_rcvs,nullptr);

	    for (int k = 0; k < N_rcvs; ++k)
	    {
                if (recv_size[k] > 0)
                {
                    auto const& cctc = thecpc.m_RcvTags->at(recv_from[k]);
                    recv_cctc[k] = &cctc;
                }
	    }

            if (actual_n_rcvs > 0) {
                stats.resize(N_rcvs);
                BL_MPI_REQUIRE( MPI_Waitall(N_rcvs, recv_reqs.dataPtr(), stats.dataPtr()) );
            }
	    
#ifdef _OPENMP
#pragma omp parallel if (FAB::isCopyOMPSafe() && thecpc.m_threadsafe_rcv)
#endif
            {
                FAB fab;
                
#ifdef _OPENMP
#pragma omp for
#endif
                for (int k = 0; k < N_rcvs; ++k)
                {
                    const char* dptr = recv_data[k];
                    if (dptr != nullptr)
                    {
                        auto const& cctc = *recv_cctc[k];
                        for (auto const& tag : cctc)
                        {
                            const Box& bx  = tag.dbox;
                            std::size_t n;
                            if (op == FabArrayBase::COPY)
                            {
                                n = (*dest)[tag.dstIndex].copyFromMem(bx,DC,NC,dptr);
                            }
                            else
                            {
                                fab.resize(bx,NC);
                                n = fab.copyFromMem(bx,0,NC,dptr);
                                (*dest)[tag.dstIndex].plus(fab,bx,bx,0,DC,NC);
                            }
                            dptr += n;
                        }
                        BL_ASSERT(dptr == recv_data[k] + recv_size[k]);
                    }
                }
	    }

#ifdef AMREX_USE_DEVICE
            Device::synchronize();
#endif
	
            if (the_recv_data) {
                amrex::The_FA_Arena()->free(the_recv_data);
            }
            else
            {
                for (auto p : recv_data) {
                    amrex::The_FA_Arena()->free(p);
                }
            }

        }
	
        if (N_snds > 0) {
            if (FabArrayBase::do_async_sends && ! thecpc.m_SndTags->empty()) {
                FabArrayBase::WaitForAsyncSends(thecpc.m_SndTags->size(),send_reqs,send_data,stats);
            }
        }

        ipass     += NC;
        SC        += NC;
        DC        += NC;
        NCompLeft -= NC;
    }

    return;
#endif /*BL_USE_MPI*/
}


template <class FAB>
template <class>
void
FabArray<FAB>::setVal (value_type val)
{
    BL_PROFILE("FabArray::setVal(val)");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter fai(*this,true); fai.isValid(); ++fai)
    {
	const Box& bx = fai.growntilebox();
	const int idx = fai.tileIndex();

        get(fai).setVal(val, bx, 0, n_comp);
    }

}

template <class FAB>
template <class>
void
FabArray<FAB>::operator= (const value_type& val)
{
    setVal(val);
}

template <class FAB>
template <class>
void
FabArray<FAB>::setVal (value_type val,
                       int        comp,
                       int        ncomp,
                       int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+ncomp <= n_comp);

    BL_PROFILE("FabArray::setVal(val,comp,ncomp,nghost)");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter fai(*this,true); fai.isValid(); ++fai)
    {
	const Box& bx = fai.growntilebox(nghost);
        get(fai).setVal(val, bx, comp, ncomp);
    }
}

template <class FAB>
template <class>
void
FabArray<FAB>::setVal (value_type val,
                       const Box& region,
                       int        comp,
                       int        ncomp,
                       int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+ncomp <= n_comp);

    BL_PROFILE("FabArray::setVal(val,region,comp,ncomp,nghost)");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter fai(*this,true); fai.isValid(); ++fai)
    {
        Box b = fai.growntilebox(nghost) & region;

        if (b.ok())
            get(fai).setVal(val, b, comp, ncomp);
    }
}

template <class FAB>
void
FabArray<FAB>::shift (const IntVect& v)
{
    clearThisBD();  // The new boxarry will have a different ID.
    for(int id(0); id < BL_SPACEDIM; ++id)
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
FabArray<FAB>::LockAllFAPointers ()
{
  typename std::map<int, std::map<int, FabArray<FAB> *> >::iterator afapIter;
  for(afapIter = FabArray<FAB>::allocatedFAPointers.begin();
      afapIter != FabArray<FAB>::allocatedFAPointers.end(); ++afapIter)
  {
    std::map<int, FabArray<FAB> *> &faPtrCachedMap = afapIter->second;
    for(typename std::map<int, FabArray<FAB> *>::iterator it = faPtrCachedMap.begin();
        it != faPtrCachedMap.end(); ++it)
    {
      FabArray<FAB> *faPtr = it->second;
      faPtr->aFAPIdLock = 1;
    }
  }
}


template <class FAB>
void
FabArray<FAB>::CheckFAPointers (bool abortOnError)
{
  int pLocked(0), pUnlocked(0);
  int allocFAPSize(FabArray<FAB>::allocatedFAPointers.size());
  int iopNumber(ParallelDescriptor::IOProcessorNumber());
  ParallelDescriptor::Bcast(&allocFAPSize, 1, iopNumber);

  if(allocFAPSize != FabArray<FAB>::allocatedFAPointers.size()) {
    if(abortOnError) {
      amrex::USleep(ParallelDescriptor::MyProcAll());
      std::cout << ParallelDescriptor::MyProcAll() << ":: **** ERROR:  allocFAPSize != "
                << "FabArray<FAB>::allocatedFAPointers.size() : " << allocFAPSize
	        << "  " << FabArray<FAB>::allocatedFAPointers.size() << std::endl;
      amrex::Abort("**** Error in FabArray<FAB>::CheckFAPointers():  bad allocFAPSize.");
    }
  }

  typename std::map<int, std::map<int, FabArray<FAB> *> >::iterator afapIter;
  for(afapIter = FabArray<FAB>::allocatedFAPointers.begin();
      afapIter != FabArray<FAB>::allocatedFAPointers.end(); ++afapIter)
  {
    std::map<int, FabArray<FAB> *> &faPtrCachedMap = afapIter->second;
    int fapCMapSize(faPtrCachedMap.size());
    ParallelDescriptor::Bcast(&fapCMapSize, 1, iopNumber);
    if(fapCMapSize != faPtrCachedMap.size()) {
      if(abortOnError) {
        amrex::USleep(ParallelDescriptor::MyProcAll());
        std::cout << ParallelDescriptor::MyProcAll() << ":: **** ERROR:  fapCMapSize != "
	          << "faPtrCachedMap.size() : " << fapCMapSize << "  "
		  << faPtrCachedMap.size() << std::endl;
        amrex::Abort("**** Error in FabArray<FAB>::CheckFAPointers():  bad fapCMapSize.");
      }
    }

    int distmapNGrids(afapIter->first-1);
    ParallelDescriptor::Bcast(&distmapNGrids, 1, iopNumber);
    if(distmapNGrids != afapIter->first-1) {
      if(abortOnError) {
        amrex::USleep(ParallelDescriptor::MyProcAll());
        std::cout << ParallelDescriptor::MyProcAll() << ":: **** ERROR:  distmapNGrids != "
	          << "afapIter->first-1 : " << distmapNGrids << "  " << afapIter->first-1 << std::endl;
        amrex::Abort("**** Error in FabArray<FAB>::CheckFAPointers():  bad distmapNGrids.");
      }
    }

    Vector<int> fapId(faPtrCachedMap.size()), dmId(faPtrCachedMap.size());
    Vector<int> nDM(faPtrCachedMap.size()), fabsAlloced(faPtrCachedMap.size());
    for(typename std::map<int, FabArray<FAB> *>::iterator it = faPtrCachedMap.begin();
        it != faPtrCachedMap.end(); ++it)
    {
      FabArray<FAB> *faPtr = it->second;
      fapId.push_back(faPtr->AllocatedFAPtrID());
      dmId.push_back(faPtr->DistributionMap().DistMapID());
      nDM.push_back(faPtr->DistributionMap().NDistMaps());
      std::set<int>::iterator nait = noallocFAPIds.find(faPtr->AllocatedFAPtrID());
      bool fabsAllocated(nait == noallocFAPIds.end());
      fabsAlloced.push_back(fabsAllocated);

      if(faPtr->aFAPIdLock) {
        ++pLocked;
      } else {
        ++pUnlocked;
      }
    }
    Vector<int> fapIdCheck, dmIdCheck, nDMCheck, fabsAllocedCheck;
    if(ParallelDescriptor::IOProcessor()) {
      fapIdCheck = fapId;
      dmIdCheck = dmId;
      nDMCheck = nDM;
      fabsAllocedCheck = fabsAlloced;
    }
    amrex::BroadcastArray(fapIdCheck, ParallelDescriptor::MyProc(), iopNumber,
                           ParallelDescriptor::Communicator());
    amrex::BroadcastArray(dmIdCheck, ParallelDescriptor::MyProc(), iopNumber,
                           ParallelDescriptor::Communicator());
    amrex::BroadcastArray(nDMCheck, ParallelDescriptor::MyProc(), iopNumber,
                           ParallelDescriptor::Communicator());
    amrex::BroadcastArray(fabsAllocedCheck, ParallelDescriptor::MyProc(), iopNumber,
                           ParallelDescriptor::Communicator());
    for(int i(0); i < fapId.size(); ++i) {
      if(fapId[i] != fapIdCheck[i]) {
        if(abortOnError) {
          std::cout << ParallelDescriptor::MyProcAll() << ":: **** ERROR:  fapId[" << i
	            << "] != fapIdCheck[" << i << "] : "
                    << fapId[i] << "  " << fapIdCheck[i] << std::endl;
          amrex::Abort("**** Error in FabArray<FAB>::CheckFAPointers():  bad fapId.");
        }
      }
    }
    for(int i(0); i < dmId.size(); ++i) {
      if(dmId[i] != dmIdCheck[i]) {
        if(abortOnError) {
          std::cout << ParallelDescriptor::MyProcAll() << ":: **** ERROR:  dmId[" << i
	            << "] != dmIdCheck[" << i << "] : "
                    << dmId[i] << "  " << dmIdCheck[i] << std::endl;
          amrex::Abort("**** Error in FabArray<FAB>::CheckFAPointers():  bad dmId.");
        }
      }
    }
    for(int i(0); i < nDM.size(); ++i) {
      if(nDM[i] != nDMCheck[i]) {
        if(abortOnError) {
          std::cout << ParallelDescriptor::MyProcAll() << ":: **** ERROR:  nDM[" << i
	            << "] != nDMCheck[" << i << "] : "
                    << nDM[i] << "  " << nDMCheck[i] << std::endl;
          amrex::Abort("**** Error in FabArray<FAB>::CheckFAPointers():  bad nDM.");
        }
      }
    }
    for(int i(0); i < fabsAlloced.size(); ++i) {
      if(fabsAlloced[i] != fabsAllocedCheck[i]) {
        if(abortOnError) {
          std::cout << ParallelDescriptor::MyProcAll() << ":: **** ERROR:  fabsAlloced[" << i
	            << "] != fabsAllocedCheck[" << i << "] : "
                    << fabsAlloced[i] << "  " << fabsAllocedCheck[i] << std::endl;
          amrex::Abort("**** Error in FabArray<FAB>::CheckFAPointers():  bad fabsAlloced.");
        }
      }
    }
  }
  int pLockedCheck(pLocked);
  ParallelDescriptor::Bcast(&pLockedCheck, 1, iopNumber);
  if(pLocked != pLockedCheck) {
    if(abortOnError) {
      std::cout << ParallelDescriptor::MyProcAll() << ":: **** ERROR:  "
                << "pLocked != pLockedCheck : " << pLocked << "  " << pLockedCheck << std::endl;
      amrex::Abort("**** Error in FabArray<FAB>::CheckFAPointers():  bad pLocked.");
    }
  }
}


template <class FAB>
void
FabArray<FAB>::PrintFAPointers ()
{
  std::cout << "|||| _in PrintFAPointers:  allocatedFAPointers.size() = "
            << FabArray<FAB>::allocatedFAPointers.size() << std::endl;
  typename std::map<int, std::map<int, FabArray<FAB> *> >::iterator afapIter;
  for(afapIter = FabArray<FAB>::allocatedFAPointers.begin();
      afapIter != FabArray<FAB>::allocatedFAPointers.end(); ++afapIter)
  {
    std::map<int, FabArray<FAB> *> &faPtrCachedMap = afapIter->second;
    std::cout << "|||| |||| _in PrintFAPointers:  faPtrCachedMap.size() distmapngrids = "
              << faPtrCachedMap.size() << "  " << afapIter->first-1 << std::endl;
    for(typename std::map<int, FabArray<FAB> *>::iterator it = faPtrCachedMap.begin();
        it != faPtrCachedMap.end(); ++it)
    {
      FabArray<FAB> *faPtr = it->second;
      std::set<int>::iterator nait = noallocFAPIds.find(faPtr->AllocatedFAPtrID());
      bool fabsAllocated(nait == noallocFAPIds.end());
      std::cout << "|||| |||| |||| _in PrintFAPointers:  FAPtrId  idLock dmID nDM fabsAlloced = " 
                << faPtr->AllocatedFAPtrID() << "  " << faPtr->aFAPIdLock
		<< "  " << faPtr->DistributionMap().DistMapID()
		<< "  " << faPtr->DistributionMap().NDistMaps()
		<< "  " << fabsAllocated << std::endl;
    }
  }
  std::cout << "|||| |||| |||| _out PrintFAPointers: ================ " << std::endl;
}


template <class FAB>
void
FabArray<FAB>::AddProcsToComp (int ioProcNumSCS, int ioProcNumAll,
                               int scsMyId, MPI_Comm scsComm)
{
  flushTileArrayCache();

  // ---- BoxArrays
  amrex::BroadcastBoxArray(boxarray, scsMyId, ioProcNumSCS, scsComm);

  // ---- DistributionMapping
  amrex::BroadcastDistributionMapping(distributionMap, scsMyId, ioProcNumSCS, scsComm, true);

  // ---- ints
  ParallelDescriptor::Bcast(&n_grow, 1, ioProcNumSCS, scsComm);
  ParallelDescriptor::Bcast(&n_comp, 1, ioProcNumSCS, scsComm);
  ParallelDescriptor::Bcast(&aFAPId, 1, ioProcNumSCS, scsComm);
  ParallelDescriptor::Bcast(&nFabArrays, 1, ioProcNumSCS, scsComm);

  // ---- bools
  bool bii(IsInitialized());
  int bInit(bii), bDAS(do_async_sends);
  ParallelDescriptor::Bcast(&bInit, 1, ioProcNumSCS, scsComm);
  ParallelDescriptor::Bcast(&bDAS, 1, ioProcNumSCS, scsComm);
  if(scsMyId != ioProcNumSCS) {
    SetInitialized(bInit);
    do_async_sends = bDAS;
  }

  // ---- IntVects
  Vector<int> allIntVects;
  if(scsMyId == ioProcNumSCS) {
    for(int i(0); i < BL_SPACEDIM; ++i)    { allIntVects.push_back(mfiter_tile_size[i]); }
    for(int i(0); i < BL_SPACEDIM; ++i)    { allIntVects.push_back(comm_tile_size[i]); }
  }
  amrex::BroadcastArray(allIntVects, scsMyId, ioProcNumSCS, scsComm);
  if(scsMyId != ioProcNumSCS) {
    int count(0);
    for(int i(0); i < BL_SPACEDIM; ++i)    { mfiter_tile_size[i] = allIntVects[count++]; }
    for(int i(0); i < BL_SPACEDIM; ++i)    { comm_tile_size[i] = allIntVects[count++]; }
  }

  // ---- add to allocatedFAPointers
  int isAllocated(0);
  if(scsMyId == ioProcNumSCS) {
    typename std::map<int, std::map<int, FabArray<FAB> *> >::iterator afapIter =
                 FabArray<FAB>::allocatedFAPointers.find(distributionMap.size());
    if(afapIter == FabArray<FAB>::allocatedFAPointers.end()) {  // ---- not found
      isAllocated = 0;
    } else {  // ---- now look for the fapid
      std::map<int, FabArray<FAB> *> &faPtrCachedMap = afapIter->second;
      if(faPtrCachedMap.find(aFAPId) != faPtrCachedMap.end()) {
        isAllocated = 1;
      }
    }
  }
  ParallelDescriptor::Bcast(&isAllocated, 1, ioProcNumSCS, scsComm);
  if(scsMyId != ioProcNumSCS) {
    if(isAllocated) {  // ---- insert into the maps
      typename std::map<int, std::map<int, FabArray<FAB> *> >::iterator afapIter =
                   FabArray<FAB>::allocatedFAPointers.find(distributionMap.size());
      if(afapIter == FabArray<FAB>::allocatedFAPointers.end()) {
        int dmapSize(distributionMap.size());
        std::map<int, FabArray<FAB> *> tempMap;
        tempMap.insert(std::make_pair(aFAPId, this));
        FabArray<FAB>::allocatedFAPointers.insert(std::make_pair(dmapSize, tempMap));

      } else {
        std::map<int, FabArray<FAB> *> &faPtrCachedMap = afapIter->second;
        faPtrCachedMap.insert(std::make_pair(aFAPId, this));
      }
    }
  }

  // ---- noallocFAPIds
  Vector<int> setArray;
  if(scsMyId == ioProcNumSCS) {
    std::set<int>::iterator it;
    for(it = noallocFAPIds.begin(); it != noallocFAPIds.end(); ++it) {
      setArray.push_back(*it);
    }
  }
  amrex::BroadcastArray(setArray, scsMyId, ioProcNumSCS, scsComm);
  if(scsMyId != ioProcNumSCS) {
      noallocFAPIds.insert(std::begin(setArray), std::end(setArray));
  }

  // ---- m_bdkey
  if(scsMyId != ioProcNumSCS) {
    m_bdkey = getBDKey();
  }

  // ---- static TACache     m_TheTileArrayCache;     // ---- cleared with flushTileArrayCache();
  // ---- static CacheStats  m_TAC_stats;             // ---- leave as is for now
  // ---- m_BD_count                                  // ---- leave as is for now
  // ---- m_FA_stats                                  // ---- leave as is for now


  // ---- unlock
  aFAPIdLock = 0;
}


// Move and replace distribution maps in all FabArrays that 
// have the same size as newDistMap.
// Note: Check for and skip over case of newDistMap = currentMap?
template <class FAB>
void
FabArray<FAB>::MoveAllFabs (const DistributionMapping &newDistMap)
{
  const Array<int> &newDistMapArray = newDistMap.ProcessorMap();
if (ParallelDescriptor::IOProcessor()) {
  std::cout << "####### start MoveAllFabs ########" << std::endl;
  std::cout << "new dist map: RefId = " << newDistMap.getRefID() << ", size = " << newDistMap.size() << std::endl;
/*
  std::cout << "####### start MoveAllFabs ##### new dm ###" << std::endl;
  for(int idm(0); idm < newDistMapArray.size(); ++idm) {
    std::cout << "newDistMapArray[" << idm << "] = " << newDistMapArray[idm] << std::endl;
  }
}*/
  std::cout << "##########################################" << std::endl;
}
  //DistributionMapping newDistMap(newDistMapArray);

//  std::map<BDKey, FabArray<FAB> *> replaceTheseDistMaps;
  if(ParallelDescriptor::IOProcessor()) {
    //std::cout << "000-----------------------------------------------------------" << std::endl;
    std::cout << "FabArray<FAB>::MoveAllFabs:  " << allocatedFAPointers.size()
              << " cached DistributionMap(s)." << std::endl;
    for(auto afapIter = allocatedFAPointers.begin(); afapIter != allocatedFAPointers.end(); ++afapIter) {
      std::cout << "afapIter->first (ngrids) = " << afapIter->first << "  :: ===================" << std::endl;
      std::map<int, FabArray<FAB> *> &afapMap = afapIter->second;
      for(auto afapMapIter = afapMap.begin(); afapMapIter != afapMap.end(); ++afapMapIter) {
	auto faptr = afapMapIter->second;
        std::cout << "afapMapIter->first (aFAPId) = " << afapMapIter->first << "  faptr = " << faptr
	          << ", BDKey = " << faptr->getBDKey() << std::endl;
      }
    }
    //std::cout << "000-----------------------------------------------------------" << std::endl;
  }

/*
  ParallelDescriptor::Barrier();

  if(ParallelDescriptor::MyProc() == 6) {
    std::cout << "666-----------------------------------------------------------" << std::endl;
    std::cout << "FabArray<FAB>::MoveAllFabs:  " << allocatedFAPointers.size()
              << " cached DistributionMap(s)." << std::endl;
    for(auto afapIter = allocatedFAPointers.begin(); afapIter != allocatedFAPointers.end(); ++afapIter) {
      std::cout << "afapIter->first (ngrids) = " << afapIter->first << "  :: ===================" << std::endl;
      std::map<int, FabArray<FAB> *> &afapMap = afapIter->second;
      for(auto afapMapIter = afapMap.begin(); afapMapIter != afapMap.end(); ++afapMapIter) {
	auto faptr = afapMapIter->second;
        std::cout << "afapMapIter->first (aFAPId) = " << afapMapIter->first << "  faptr = " << faptr
	          << ", BDKey = " << faptr->getBDKey() << std::endl;
      }
    }
    std::cout << "666-----------------------------------------------------------" << std::endl;
  }
*/
  // Get FabArrays with the appropriate size (that matches the new map size).
  typename std::map<int, std::map<int, FabArray<FAB> *> >::iterator afapIter = 
               FabArray<FAB>::allocatedFAPointers.find(newDistMapArray.size());
  if(afapIter == FabArray<FAB>::allocatedFAPointers.end()) {
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "FabArray<FAB>::MoveAllFabs:  no allocated pointers for new "
                << "distribution map array.\n";
    }
  } else {
    std::map<int, FabArray<FAB> *> &faPtrCachedMap = afapIter->second;
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "............................................." << std::endl
                << "FabArray<FAB>::MoveAllFabs:  moving FABs for "
                << faPtrCachedMap.size() << " FabArray(s).\n";
    }
/*
    int cached_size = faPtrCachedMap.size();
    int size_min = cached_size;
    int size_max = cached_size;
    ParallelDescriptor::ReduceIntMin(size_min); 
    ParallelDescriptor::ReduceIntMax(size_max);    
    std::cout << ParallelDescriptor::MyProc() << " : Size of faPtrCachedMap: " << cached_size << std::endl;
    BL_ASSERT(size_min == size_max); 
*/
    int loop_num = 0;
    for(typename std::map<int, FabArray<FAB> *>::iterator it = faPtrCachedMap.begin();
        it != faPtrCachedMap.end(); ++it)
    {
      loop_num++;
      FabArray<FAB> *faPtr = it->second;
//      std::cout << ParallelDescriptor::MyProc() << ", " << loop_num << " : Entering faPtr->ok() " << std::endl;
      if(faPtr->ok() == false) {
        amrex::Abort("**** Error 0 in MoveAllFabs:  faPtr not ok");
      }
//      std::cout << ParallelDescriptor::MyProc() << ", " << loop_num << " : Exiting faPtr->ok() " << std::endl;

      faPtr->flushFPinfo();

      // Move Fabs in Array based on new Distribution Map
      int nFabsMoved = faPtr->MoveFabs(newDistMapArray);

      // Update distribution map and ID key for FabArray
      amrex::Print() << "RRRRRRRR::  replacing:  " << faPtr->DistributionMap().getRefID() << " with : " << newDistMap.getRefID() << std::endl;
      //const DistributionMapping::Ref *refidPtr = faPtr->DistributionMap().getRefID().dataPtr();
      //faPtr->ModifyDistributionMap().ReplaceCachedProcessorMap(newDistMapArray);
      faPtr->ModifyDistributionMap() = newDistMap;
      faPtr->updateBDKey();

      //replaceTheseDistMaps[faPtr->getBDKey()] = faPtr;

      // REMOVE! Requires Sync!! 
      ParallelDescriptor::ReduceIntSum(nFabsMoved);
      amrex::Print() << ParallelDescriptor::MyProc() << " :FabArray<FAB>::MoveAllFabs:  moved " << nFabsMoved
                     << " FAB(s) for faPtr = " << faPtr << '\n';
    }

    // ---- flush caches
    flushCPCache();
    flushFBCache();
    flushTileArrayCache();

    for(typename std::map<int, FabArray<FAB> *>::iterator it = faPtrCachedMap.begin();
        it != faPtrCachedMap.end(); ++it)
    {
      FabArray<FAB> *faPtr = it->second;
      faPtr->flushFPinfo();
      faPtr->flushCFinfo();
    }

//    for(auto rTDMI = replaceTheseDistMaps.begin(); rTDMI != replaceTheseDistMaps.end(); ++rTDMI) {
//      FabArray<FAB> *faPtr = rTDMI->second;
//      amrex::Print() << "RRRRRRRR::  replacing:  " << faPtr->DistributionMap() << std::endl;
//      const DistributionMapping::Ref *refidPtr = faPtr->DistributionMap().getRefID().dataPtr();
//      amrex::Print() << "DDDDDDDD::  refidPtr =  " << refidPtr << std::endl;
//      //faPtr->ModifyDistributionMap().ReplaceCachedProcessorMap(newDistMapArray);
//      faPtr->ModifyDistributionMap() = newDistMap;
//      faPtr->updateBDKey();
//    }

    // ---- check the moved fabarrays
    for(typename std::map<int, FabArray<FAB> *>::iterator it = faPtrCachedMap.begin();
        it != faPtrCachedMap.end(); ++it)
    {
      FabArray<FAB> *faPtr = it->second;
      if(faPtr->ok() == false) {
        amrex::Abort("**** Error 2 in MoveAllFabs:  faPtr not ok");
      }
    }

  }

  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "FabArray<FAB>::MoveAllFabs:  noallocFAPIds.size() = " << noallocFAPIds.size() << std::endl;
    std::set<int>::iterator nait;
    for(nait = noallocFAPIds.begin(); nait != noallocFAPIds.end(); ++nait) {
        std::cout << "FabArray<FAB>::MoveAllFabs:  noallocFAPIds = " << *nait << std::endl;
    }

    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "FabArray<FAB>::MoveAllFabs:  " << allocatedFAPointers.size()
              << " cached DistributionMap(s)." << std::endl;
    for(auto it = allocatedFAPointers.begin(); it != allocatedFAPointers.end(); ++it) {
      std::cout << "it->first (ngrids) = " << it->first << "  :: ===================" << std::endl;
      std::map<int, FabArray<FAB> *> &afapMap = it->second;
      for(auto afapMapIter = afapMap.begin(); afapMapIter != afapMap.end(); ++afapMapIter) {
	auto fap = afapMapIter->second;
        std::cout << "afapMapIter->first (aFAPId) = " << afapMapIter->first << "  faptr = " << afapMapIter->second 
	          << ", BDKey = " << fap->getBDKey() << std::endl;
      }
    }
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "######### end MoveAllFabs ################" << std::endl;
    }

}


// Replace all distribution maps in FabArrays that are currently using
// oldDistMap with newDistMap and move Fabs approriately.
template <class FAB>
void
//FabArray<FAB>::MoveAllFabs (const Vector<int> &newDistMapArray)
FabArray<FAB>::MoveAllFabs (const DistributionMapping &oldDistMap, const DistributionMapping &newDistMap)
{

  // If the maps are identical, skip.
  if (oldDistMap == newDistMap)
  { return; }

  int myProc(ParallelDescriptor::MyProc());
  const Vector<int> &newDistMapArray = newDistMap.ProcessorMap();
if(myProc == 0) {
  std::cout << "####### start MoveAllFabs ##### new dm ###" << std::endl;
  for(int idm(0); idm < newDistMapArray.size(); ++idm) {
    std::cout << "newDistMapArray[" << idm << "] = " << newDistMapArray[idm] << std::endl;
  }
  std::cout << "##########################################" << std::endl;
}

  std::map<BDKey, FabArray<FAB> *> replaceTheseDistMaps;
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "FabArray<FAB>::MoveAllFabs:  " << allocatedFAPointers.size()
              << " cached DistributionMap(s)." << std::endl;
    for(auto afapIter = allocatedFAPointers.begin(); afapIter != allocatedFAPointers.end(); ++afapIter) {
      std::cout << "afapIter->first (ngrids) = " << afapIter->first << "  :: ===================" << std::endl;
      std::map<int, FabArray<FAB> *> &afapMap = afapIter->second;
      for(auto afapMapIter = afapMap.begin(); afapMapIter != afapMap.end(); ++afapMapIter) {
	auto faptr = afapMapIter->second;
        std::cout << "afapMapIter->first (aFAPId) = " << afapMapIter->first << "  faptr = " << faptr << std::endl;
	std::cout << "faptr->DistributionMap() = " << faptr->DistributionMap() << std::endl;
      }
    }
    std::cout << "-----------------------------------------------------------" << std::endl;
  }

  // Select subset of allocatedFabArrays that are the correct size for the distribution map.
  typename std::map<int, std::map<int, FabArray<FAB> *> >::iterator afapIter = 
               FabArray<FAB>::allocatedFAPointers.find(newDistMapArray.size());
  if(afapIter == FabArray<FAB>::allocatedFAPointers.end()) {
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "FabArray<FAB>::MoveAllFabs:  no allocated pointers for new "
                << "distribution map array.\n";
    }
  } else {
    std::map<int, FabArray<FAB> *> &faPtrCachedMap = afapIter->second;
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "FabArray<FAB>::MoveAllFabs:  moving FABs for "
                << faPtrCachedMap.size() << " FabArray(s).\n";
    }

    for(typename std::map<int, FabArray<FAB> *>::iterator it = faPtrCachedMap.begin();
        it != faPtrCachedMap.end(); ++it)
    {
      FabArray<FAB> *faPtr = it->second;
      if(faPtr->ok() == false) {
        amrex::Abort("**** Error 0 in MoveAllFabs:  faPtr not ok");
      }

      // Find the appropriate FabArrays that use the old Map.
      if ( faPtr->ModifyDistributionMap() == oldDistMap )
      {

        faPtr->flushFPinfo();
        // Move the fabs based on the new map.
        int nFabsMoved = faPtr->MoveFabs(newDistMapArray);

        // Update the distribution map and the unique fab array ID.
        amrex::Print() << "RRRRRRRR::  replacing:  " << faPtr->DistributionMap() << std::endl;
        const DistributionMapping::Ref *refidPtr = faPtr->DistributionMap().getRefID().dataPtr();
        amrex::Print() << "DDDDDDDD::  refidPtr =  " << refidPtr << std::endl;
        faPtr->ModifyDistributionMap() = newDistMap;
        faPtr->updateBDKey();

        // REMOVE! Requires Sync!! 
        ParallelDescriptor::ReduceIntSum(nFabsMoved);
        amrex::Print() << "FabArray<FAB>::MoveAllFabs:  moved " << nFabsMoved
                       << " FAB(s) for faPtr = " << faPtr << '\n';
      }
    }



    // ---- flush caches
    flushCPCache();
    flushFBCache();
    flushTileArrayCache();

    for(typename std::map<int, FabArray<FAB> *>::iterator it = faPtrCachedMap.begin();
        it != faPtrCachedMap.end(); ++it)
    {
      FabArray<FAB> *faPtr = it->second;
      faPtr->flushFPinfo();
      faPtr->flushCFinfo();
    }

    // ---- check the moved fabarrays
    for(typename std::map<int, FabArray<FAB> *>::iterator it = faPtrCachedMap.begin();
        it != faPtrCachedMap.end(); ++it)
    {
      FabArray<FAB> *faPtr = it->second;
      if(faPtr->ok() == false) {
        amrex::Abort("**** Error 2 in MoveAllFabs:  faPtr not ok");
      }
    }

  }

  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "FabArray<FAB>::MoveAllFabs:  noallocFAPIds.size() = " << noallocFAPIds.size() << std::endl;
    std::set<int>::iterator nait;
    for(nait = noallocFAPIds.begin(); nait != noallocFAPIds.end(); ++nait) {
        std::cout << "FabArray<FAB>::MoveAllFabs:  noallocFAPIds = " << *nait << std::endl;
    }
  }
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "FabArray<FAB>::MoveAllFabs:  " << allocatedFAPointers.size()
              << " cached DistributionMap(s)." << std::endl;
    for(auto it = allocatedFAPointers.begin(); it != allocatedFAPointers.end(); ++it) {
      std::cout << "it->first (ngrids) = " << it->first << "  :: ===================" << std::endl;
      std::map<int, FabArray<FAB> *> &afapMap = it->second;
      for(auto afapMapIter = afapMap.begin(); afapMapIter != afapMap.end(); ++afapMapIter) {
	auto fap = afapMapIter->second;
        std::cout << "afapMapIter->first (aFAPId) = " << afapMapIter->first << "  faptr = " << afapMapIter->second << std::endl; 
	std::cout << "fap->DistributionMap() = " << fap->DistributionMap() << std::endl;
      }
    }
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  }
if(myProc == 0) {
  std::cout << "######### end MoveAllFabs ################" << std::endl;
}
}

template <class FAB>
int
FabArray<FAB>::MoveFabs (const Vector<int> &newDistMapArray)
{
#if BL_USE_MPI
  BL_PROFILE("FabArray<FAB>::MoveFabs()");

  int myProc(ParallelDescriptor::MyProc());

  // ---- check validity of newDistMapArray
  if(newDistMapArray.size() != distributionMap.size()) {
    std::cout << "ndma.size  dm.size = " << newDistMapArray.size() << "  "
              << distributionMap.size() << std::endl;
    amrex::Abort("**** Error:  bad newDistMap:0");
  }
  for(int idm(0); idm < newDistMapArray.size(); ++idm) {
    if(newDistMapArray[idm] < 0 ||
       newDistMapArray[idm] > ParallelDescriptor::NProcs() - 1)
    {
      std::cout << "**** idm ndm[idm] np = " << newDistMapArray[idm] << "  "
                << ParallelDescriptor::NProcs() << std::endl;
      amrex::Abort("**** Error:  bad newDistMap:2");
    }
  }
  if(newDistMapArray == distributionMap.ProcessorMap()) {
    return 0;
  }

  // ---- try to resolve fabarrays that were created without allocating memory
  std::set<int>::iterator it = noallocFAPIds.find(aFAPId);
  bool fabsNotAllocated(it != noallocFAPIds.end());
  if(fabsNotAllocated) {
    if(ok()) {  // ---- this fabarray was created without memory allocated then allocated later
      noallocFAPIds.erase(aFAPId);
    } else {
      return 0;
    }
  } else {
    BL_ASSERT(ok());
  }

  // ---- determine which fabs to move
  std::vector<FABMoves> fabMoves;
  for(int iM(0); iM < distributionMap.size(); ++iM) {
    if(newDistMapArray[iM] != distributionMap[iM]) {
      FABMoves moveThisFab;
      moveThisFab.distMapIndex = iM;
      moveThisFab.fromRank     = distributionMap[iM];
      moveThisFab.toRank       = newDistMapArray[iM];
      moveThisFab.seqNum       = ParallelDescriptor::SeqNum();
      fabMoves.push_back(moveThisFab);
    }
  }

  // -- save the original index values to preserve order
  std::map<int, FAB *> tempIndexFABs;
  typename std::map<int, FAB *>::iterator tIFiter;
  if(fabMoves.size() > 0) {  // ---- there are fabs on this proc to send || recv | both
    if(indexArray.size() != m_fabs_v.size()) {
      amrex::Abort("**** Error:  indexArray size != m_fabs_v.size()");
    }
    for(int iim(0); iim < indexArray.size(); ++iim) {
      tIFiter = tempIndexFABs.find(indexArray[iim]);
      if(tIFiter == tempIndexFABs.end()) {
	tempIndexFABs.insert(std::pair<int, FAB *>(indexArray[iim], m_fabs_v[iim]));
      } else {
        amrex::Abort("**** Error:  index already in indexArray.");
      }
    }
  }

  // ---- move the fabs
  Vector<MPI_Request> recvReqs, sendReqs;
  int nFabsSent(0);

  for(int imoves(0); imoves < fabMoves.size(); ++imoves) {
    FABMoves &moveThisFab = fabMoves[imoves];
    int dmi(moveThisFab.distMapIndex);

    if(myProc == moveThisFab.toRank) {   // ---- receive fab(s)
      const Box &tmpbox = fabbox(dmi);
      FAB *fabPtr = new FAB(tmpbox, n_comp);
      tempIndexFABs[dmi] = fabPtr;  // ---- add to map

      recvReqs.push_back(ParallelDescriptor::Arecv(fabPtr->dataPtr(),
                                                   tmpbox.numPts() * n_comp,
                                                   moveThisFab.fromRank,
						   moveThisFab.seqNum).req());
    }

    if(myProc == moveThisFab.fromRank) {    // ---- send fab(s)
      ++nFabsSent;
      FAB *fabPtr;
      tIFiter = tempIndexFABs.find(dmi);
      if(tIFiter == tempIndexFABs.end()) {
	std::cout << myProc << ":: fromRank this = " << moveThisFab.fromRank << "  " << this << std::endl;
        amrex::Abort("**** Error:  index not in tempIndexFABs.");
      } else {
        fabPtr = tIFiter->second;
      }
      BL_ASSERT(fabPtr->nComp() == n_comp);

      if(FabArrayBase::do_async_sends) {
        sendReqs.push_back(ParallelDescriptor::Asend(fabPtr->dataPtr(),
                                                     fabPtr->box().numPts() * n_comp,
		                                     moveThisFab.toRank,
		                                     moveThisFab.seqNum).req());
      } else {
        ParallelDescriptor::Send(fabPtr->dataPtr(),
                                 fabPtr->box().numPts() * n_comp,
		                 moveThisFab.toRank,
		                 moveThisFab.seqNum);
      }
    }
  }

  // ---- wait for all the data to move
  // ---- we could defer Waitall for multiple calls to MoveFabs with a multi-step process
  Vector<MPI_Status>  recvStats(recvReqs.size()), sendStats(sendReqs.size());

  if(recvReqs.size() > 0) {
    BL_MPI_REQUIRE( MPI_Waitall(recvReqs.size(), recvReqs.dataPtr(), recvStats.dataPtr()) );
  }

  if(FabArrayBase::do_async_sends) {
    if(sendReqs.size() > 0) {
      BL_MPI_REQUIRE( MPI_Waitall(sendReqs.size(), sendReqs.dataPtr(), sendStats.dataPtr()) );
    }
  }

  // ---- delete moved fabs and data from the temporary map
  for(int imoves(0); imoves < fabMoves.size(); ++imoves) {
    FABMoves &moveThisFab = fabMoves[imoves];
    if(myProc == moveThisFab.fromRank) {  // ---- delete sent fab
      tIFiter = tempIndexFABs.find(moveThisFab.distMapIndex);
      if(tIFiter == tempIndexFABs.end()) {
        amrex::Abort("**** Error:  index not in tempIndexFABs when deleting.");
      } else {
	delete tIFiter->second;
	tempIndexFABs.erase(tIFiter);
      }
    }
  }

  // ---- reconstruct the index and fab vectors
  indexArray.clear();
  ownership.clear();
  m_fabs_v.clear();
  for(tIFiter = tempIndexFABs.begin(); tIFiter != tempIndexFABs.end(); ++tIFiter) {
    indexArray.push_back(tIFiter->first);
    ownership.push_back(myProc == newDistMapArray[tIFiter->first]);
    m_fabs_v.push_back(tIFiter->second);
  }

  amrex::Print() << "(amrex::FabArray::MoveFabs):  updating bdkey for faPtr =  " << this << std::endl;
  updateBDKey();

  return nFabsSent;
#else
  return 0;
#endif
}

template <class FAB>
void
FabArray<FAB>::FillBoundary (bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    if ( n_grow > 0 ) {
	FillBoundary_nowait(0, nComp(), Periodicity::NonPeriodic(), cross);
	FillBoundary_finish();
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundary (const Periodicity& period, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    if ( n_grow > 0 ) {
	FillBoundary_nowait(0, nComp(), period, cross);
	FillBoundary_finish();
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundary (int scomp, int ncomp, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    if ( n_grow > 0 ) {
	FillBoundary_nowait(scomp, ncomp, Periodicity::NonPeriodic(), cross);
	FillBoundary_finish();
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundary (int scomp, int ncomp, const Periodicity& period, bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");
    if ( n_grow > 0 ) {
	FillBoundary_nowait(scomp, ncomp, period, cross);
	FillBoundary_finish();
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundary_nowait (bool cross)
{
    FillBoundary_nowait(0, nComp(), Periodicity::NonPeriodic(), cross);
}

template <class FAB>
void
FabArray<FAB>::FillBoundary_nowait (const Periodicity& period, bool cross)
{
    FillBoundary_nowait(0, nComp(), period, cross);
}

template <class FAB>
void
FabArray<FAB>::FillBoundary_nowait (int scomp, int ncomp, bool cross)
{
    FillBoundary_nowait(scomp, ncomp, Periodicity::NonPeriodic(), cross);
}

template <class FAB>
void
FabArray<FAB>::EnforcePeriodicity (const Periodicity& period)
{
    BL_PROFILE("FabArray::EnforcePeriodicity");
    if (period.isAnyPeriodic()) {
	FBEP_nowait(0, nComp(), period, false, true);
	FillBoundary_finish(); // unsafe unless isAnyPeriodic()
    }
}

template <class FAB>
void
FabArray<FAB>::EnforcePeriodicity (int scomp, int ncomp, const Periodicity& period)
{
    BL_PROFILE("FabArray::EnforcePeriodicity");
    if (period.isAnyPeriodic()) {
	FBEP_nowait(scomp, ncomp, period, false, true);
	FillBoundary_finish(); // unsafe unless isAnyPeriodic()
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundary_nowait (int scomp, int ncomp, const Periodicity& period, bool cross)
{
    FBEP_nowait(scomp, ncomp, period, cross);
}

template <class FAB>
void
FabArray<FAB>::FBEP_nowait (int scomp, int ncomp, const Periodicity& period, bool cross,
			    bool enforce_periodicity_only)
{
    fb_cross = cross;
    fb_epo   = enforce_periodicity_only;
    fb_scomp = scomp;
    fb_ncomp = ncomp;
    fb_period = period;

    bool work_to_do;
    if (enforce_periodicity_only) {
	work_to_do = period.isAnyPeriodic();
    } else {
	work_to_do = n_grow > 0;
    }
    if (!work_to_do) return;

    const FB& TheFB = getFB(period, cross, enforce_periodicity_only);

    if (ParallelDescriptor::NProcs() == 1)
    {
        //
        // There can only be local work to do.
        //
	int N_loc = (*TheFB.m_LocTags).size();
#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe() && TheFB.m_threadsafe_loc)
#endif
	for (int i=0; i<N_loc; ++i)
        {
            const CopyComTag& tag = (*TheFB.m_LocTags)[i];

            BL_ASSERT(distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc());
            BL_ASSERT(distributionMap[tag.srcIndex] == ParallelDescriptor::MyProc());

            get(tag.dstIndex).copy(get(tag.srcIndex),tag.sbox,scomp,tag.dbox,scomp,ncomp);
        }

        return;
    }
   
#ifdef BL_USE_MPI

#if defined(BL_USE_UPCXX)
    ParallelDescriptor::Mode.set_upcxx_mode();
    ParallelDescriptor::Mode.incr_upcxx();
#endif

#if defined(BL_USE_MPI3)
    BL_ASSERT(FAB::preAllocatable());
#else
    BL_ASSERT(!ParallelDescriptor::MPIOneSided());
#endif

    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    int SeqNum, preSeqNum;
    {
	ParallelDescriptor::Color mycolor = this->color();
	if (mycolor == ParallelDescriptor::DefaultColor()) {
            if (!FAB::preAllocatable()) {
                preSeqNum = ParallelDescriptor::SeqNum();
            }
	    SeqNum = ParallelDescriptor::SeqNum();
	} else if (mycolor == ParallelDescriptor::SubCommColor()) {
            if (!FAB::preAllocatable()) {
                preSeqNum = ParallelDescriptor::SubSeqNum();
            }
	    SeqNum = ParallelDescriptor::SubSeqNum();
	}
	// else I don't have any data and my SubSeqNum() should not be called.
    }

    const int N_locs = TheFB.m_LocTags->size();
    const int N_rcvs = TheFB.m_RcvTags->size();
    const int N_snds = TheFB.m_SndTags->size();

    if (N_locs == 0 && N_rcvs == 0 && N_snds == 0)
        // No work to do.
        return;

    //
    // Before we post recv, let's preprocess sends in case FAB is not preAllocatable
    //
    Vector<char*> &                     send_data = fb_send_data;
    Vector<int>                         send_size;
    Vector<int>                         send_rank;
    Vector<MPI_Request>&                send_reqs = fb_send_reqs;
    Vector<const CopyComTagsContainer*> send_cctc;
    Vector<Vector<int> >                 indv_send_size;
    Vector<MPI_Request>                 pre_reqs;

    fb_tag = SeqNum;

#if defined BL_USE_UPCXX || defined BL_USE_MPI3 
    int actual_n_snds = 0;
#endif
    if (N_snds > 0)
    {
        fb_send_data.clear();
        fb_send_reqs.clear();

	send_data.reserve(N_snds);
	send_size.reserve(N_snds);
	send_rank.reserve(N_snds);
        send_reqs.reserve(N_snds);
	send_cctc.reserve(N_snds);
        indv_send_size.reserve(N_snds);

        for (auto const& kv : *TheFB.m_SndVols)
        {
            Vector<int> iss;                
            auto const& cctc = TheFB.m_SndTags->at(kv.first);

            std::size_t nbytes = 0;
            if (FAB::preAllocatable())
            {
                for (auto const& cct : kv.second)
                {
                    nbytes += (*this)[cct.srcIndex].nBytes(cct.sbox,scomp,ncomp);
                }
            }
            else
            {
                for (auto const& tag : cctc)
                {
                    std::size_t b = (*this)[tag.srcIndex].nBytes(tag.sbox,scomp,ncomp);
                    nbytes += b;
                    iss.push_back(static_cast<int>(b));
                }
            }
            
            BL_ASSERT(nbytes < std::numeric_limits<int>::max());
            
            char* data = nullptr;
            if (nbytes > 0)
            {
                data = static_cast<char*>
#ifdef BL_USE_UPCXX
                    (BLPgas::alloc(nbytes));
#else
                    (amrex::The_FA_Arena()->alloc(nbytes));
#endif
            }
                    
            send_data.push_back(data);
            send_size.push_back(static_cast<int>(nbytes));
            send_rank.push_back(kv.first);
            send_reqs.push_back(MPI_REQUEST_NULL);
            send_cctc.push_back(&cctc);
            indv_send_size.push_back(std::move(iss));
        }

#if defined BL_USE_UPCXX || defined BL_USE_MPI3         
        actual_n_snds = N_snds - std::count(send_size.begin(), send_size.end(), 0);
#endif
    }

    if (!FAB::preAllocatable())
    {
        pre_reqs.resize(N_snds,MPI_REQUEST_NULL);
        for (int j=0; j<N_snds; ++j)
        {
            pre_reqs[j] = ParallelDescriptor::Asend(indv_send_size[j].data(),
                                                    indv_send_size[j].size(),
                                                    send_rank[j],preSeqNum).req();
        }
    }

    //
    // Post rcvs. Allocate one chunk of space to hold'm all.
    //
#ifdef BL_USE_MPI3
    MPI_Group tgroup, rgroup, sgroup;
    if (ParallelDescriptor::MPIOneSided()) {
	MPI_Comm_group(ParallelDescriptor::Communicator(), &tgroup);
    }
#endif

    fb_the_recv_data = nullptr;

    if (N_rcvs > 0) {
#ifdef BL_USE_UPCXX
	PostRcvs_PGAS(*TheFB.m_RcvVols, fb_the_recv_data, fb_recv_data,
                      fb_recv_size, fb_recv_from,
                      scomp, ncomp, SeqNum, &BLPgas::fb_recv_event);
#else
	if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
	    PostRcvs_MPI_Onesided(*TheFB.m_RcvVols, fb_the_recv_data, fb_recv_data,
                                  fb_recv_size, fb_recv_from, fb_recv_reqs, fb_recv_disp,
                                  scomp, ncomp, SeqNum, ParallelDescriptor::fb_win);
	    MPI_Group_incl(tgroup, fb_recv_from.size(), fb_recv_from.dataPtr(), &rgroup);
	    MPI_Win_post(rgroup, 0, ParallelDescriptor::fb_win);
#endif
	} else {
	    PostRcvs(*TheFB.m_RcvVols, *TheFB.m_RcvTags,
                     fb_recv_data, fb_recv_size, fb_recv_from, fb_recv_reqs,
                     scomp, ncomp, SeqNum, preSeqNum);
	}
#endif
    }

    //
    // Post send's
    //
    if (N_snds > 0)
    {
#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe())
#endif
	for (int j=0; j<N_snds; ++j)
	{
            char* dptr = send_data[j];
            if (dptr != nullptr)
            {
                auto const& cctc = *send_cctc[j];
                for (auto const& tag : cctc)
                {
                    const Box& bx = tag.sbox;
                    auto n = (*this)[tag.srcIndex].copyToMem(bx,scomp,ncomp,dptr);
                    dptr += n;
                }
                BL_ASSERT(dptr == send_data[j] + send_size[j]);
            }
	}

#ifdef AMREX_USE_DEVICE
        Device::synchronize();
#endif

#ifdef BL_USE_UPCXX

	BLPgas::fb_send_counter = 0;

	for (int i=0; i<N_snds; ++i) {
            if (send_size[i] > 0) {
                BLPgas::Send(upcxx::global_ptr<void>((void *)send_data[i], upcxx::myrank()),
                             send_rank[i], send_size[i], SeqNum,
                             &BLPgas::fb_send_event, &BLPgas::fb_send_counter);
            }
        }

	// Need to make sure at least half of the sends have been started
	while (BLPgas::fb_send_counter < actual_n_snds) {
	    upcxx::advance();
        }

#else  // MPI

	if (ParallelDescriptor::MPIOneSided())
	{
#if defined(BL_USE_MPI3)
	    Vector<MPI_Aint> send_disp(N_snds,0);
            send_reqs.clear();
	    
	    for (int i=0; i<N_snds; ++i) {
                if (send_size[i] > 0) {
                    send_reqs.push_back(ParallelDescriptor::Arecv
                                        (&send_disp[i],1,send_rank[i],SeqNum).req());
                }
	    }
		
	    MPI_Group_incl(tgroup, N_snds, send_rank.dataPtr(), &sgroup);
	    MPI_Win_start(sgroup,0,ParallelDescriptor::fb_win);
	    
	    int send_counter = 0;	
	    while (send_counter < actual_n_snds) {
		MPI_Status status;
		int index;
		
		MPI_Waitany(N_snds, send_reqs.dataPtr(), &index, &status);
		
		BL_ASSERT(status.MPI_TAG == SeqNum);
		BL_ASSERT(status.MPI_SOURCE == send_rank[index]);

                MPI_Put(send_data[index], send_size[index], MPI_CHAR, send_rank[index],
                        send_disp[index], send_size[index], MPI_CHAR, ParallelDescriptor::fb_win);
		
		++send_counter;
	    }
#endif // BL_USE_MPI3
	} 
	else 
	{
            int send_counter = 0;
            while (send_counter < N_snds)
            {
                int j;
                
                if (FAB::preAllocatable())
                {
                    j = send_counter;
                }
                else
                {
                    MPI_Status status;
                    MPI_Waitany(N_snds, pre_reqs.data(), &j, &status);
                }
                        
                BL_ASSERT(send_size[j] > 0);

                send_reqs[j] = ParallelDescriptor::Asend
                    (send_data[j],send_size[j],send_rank[j],SeqNum).req();

                ++send_counter;
            }
	}
#endif // #ifdef BL_USE_UPCXX #else 
    }

    if (ParallelDescriptor::MPIOneSided()) {
#ifdef BL_USE_MPI3
	MPI_Group_free(&tgroup);
	if (N_rcvs > 0) MPI_Group_free(&rgroup);
	if (N_snds > 0) MPI_Group_free(&sgroup);
#endif
    }

    //
    // Do the local work.  Hope for a bit of communication/computation overlap.
    //
    if (ParallelDescriptor::TeamSize() > 1 && TheFB.m_threadsafe_loc)
    {
#ifdef BL_USE_TEAM
#ifdef _OPENMP
#pragma omp parallel if (FAB::isCopyOMPSafe())
#endif
	ParallelDescriptor::team_for(0, N_locs, [&] (int i) 
        {
	    const auto& tag = (*TheFB.m_LocTags)[i];

	    BL_ASSERT(ParallelDescriptor::sameTeam(distributionMap[tag.dstIndex]));
	    BL_ASSERT(ParallelDescriptor::sameTeam(distributionMap[tag.srcIndex]));

	    get(tag.dstIndex).copy(get(tag.srcIndex),tag.sbox,scomp,tag.dbox,scomp,ncomp);
	});
#endif
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe() && TheFB.m_threadsafe_loc)
#endif
	for (int i=0; i<N_locs; ++i)
	{
	    const CopyComTag& tag = (*TheFB.m_LocTags)[i];

	    BL_ASSERT(ParallelDescriptor::sameTeam(distributionMap[tag.dstIndex]));
	    BL_ASSERT(ParallelDescriptor::sameTeam(distributionMap[tag.srcIndex]));
	    
	    if (distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc()) {
		get(tag.dstIndex).copy(get(tag.srcIndex),tag.sbox,scomp,tag.dbox,scomp,ncomp);
	    }
	}
    }
#endif /*BL_USE_MPI*/
}

template <class FAB>
void
FabArray<FAB>::FillBoundary_finish ()
{
    if ( n_grow <= 0 && !fb_epo ) return; // For epo (Enforce Periodicity Only), there may be no ghost cells.

    if (ParallelDescriptor::NProcs() == 1) return;

#ifdef BL_USE_MPI

#if defined(BL_USE_UPCXX)
    ParallelDescriptor::Mode.set_upcxx_mode();
    ParallelDescriptor::Mode.decr_upcxx();
#endif

#if !defined(BL_USE_MPI3)
    BL_ASSERT(!ParallelDescriptor::MPIOneSided());
#endif

    const FB& TheFB = getFB(fb_period,fb_cross,fb_epo);

    const int N_rcvs = TheFB.m_RcvTags->size();
    const int N_snds = TheFB.m_SndTags->size();

    int actual_n_rcvs = N_rcvs - std::count(fb_recv_data.begin(), fb_recv_data.end(), nullptr);

#ifdef BL_USE_UPCXX
    if (actual_n_rcvs > 0) BLPgas::fb_recv_event.wait();
#else 
    if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
	if (N_snds > 0) MPI_Win_complete(ParallelDescriptor::fb_win);
	if (N_rcvs > 0) MPI_Win_wait    (ParallelDescriptor::fb_win);
#endif
    } else {
	if (actual_n_rcvs > 0) {
	    Vector<MPI_Status> stats(N_rcvs);
	    BL_MPI_REQUIRE( MPI_Waitall(N_rcvs, fb_recv_reqs.dataPtr(), stats.dataPtr()) );
	    if (!CheckRcvStats(stats, fb_recv_size, MPI_CHAR, fb_tag))
            {
                amrex::Abort("FillBoundary_finish failed with wrong message size");
            }
	}
    }
#endif

    if (N_rcvs > 0)
    {
	Vector<const CopyComTagsContainer*> recv_cctc(N_rcvs,nullptr);

	for (int k = 0; k < N_rcvs; k++) 
	{
            if (fb_recv_size[k] > 0)
            {
                auto const& cctc = TheFB.m_RcvTags->at(fb_recv_from[k]);
                recv_cctc[k] = &cctc;
            }
	}	

#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe() && TheFB.m_threadsafe_rcv)
#endif

	for (int k = 0; k < N_rcvs; k++) 
	{
	    const char* dptr = fb_recv_data[k];
            if (dptr != nullptr)
            {
                auto const& cctc = *recv_cctc[k];
                for (auto const& tag : cctc)
                {
                    const Box& bx  = tag.dbox;
                    std::size_t n = (*this)[tag.dstIndex].copyFromMem(bx,fb_scomp,fb_ncomp,dptr);
                    dptr += n;
                }

                BL_ASSERT(dptr == fb_recv_data[k] + fb_recv_size[k]);
            }
	}

#ifdef AMREX_USE_DEVICE
        Device::synchronize();
#endif

        if (fb_the_recv_data)
        {
#ifdef BL_USE_UPCXX
            BLPgas::free(fb_the_recv_data);
#else
            if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
                MPI_Win_detach(ParallelDescriptor::fb_win, fb_the_recv_data);
#endif
            }
	    amrex::The_FA_Arena()->free(fb_the_recv_data);
#endif
	}
        else
        {
            for (auto p : fb_recv_data) {
                amrex::The_FA_Arena()->free(p);
            }
        }        
    }

    if (N_snds > 0) {
#ifdef BL_USE_UPCXX
        FabArrayBase::WaitForAsyncSends_PGAS(N_snds,fb_send_data,
					     &BLPgas::fb_send_event,
					     &BLPgas::fb_send_counter);
#else
	if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
	    for (int i = 0; i < N_snds; ++i) {
		if (fb_send_data[i]) amrex::The_FA_Arena()->free(fb_send_data[i]);
            }
#endif
        } else {
            Vector<MPI_Status> stats;
            FabArrayBase::WaitForAsyncSends(N_snds,fb_send_reqs,fb_send_data,stats);
        }
#endif
    }

#ifdef BL_USE_TEAM
    ParallelDescriptor::MyTeam().MemoryBarrier();
#endif

#endif // MPI
}

#ifdef BL_USE_UPCXX
template <class FAB>
void
FabArray<FAB>::PostRcvs_PGAS (const MapOfCopyComTagContainers&  m_RcvVols,
                              char*&                            the_recv_data,
                              Vector<char*>&                     recv_data,
                              Vector<int>&                       recv_size,
                              Vector<int>&                       recv_from,
                              int                               icomp,
                              int                               ncomp,
                              int                               SeqNum,
                              upcxx::event*                     recv_event)
{
    recv_data.clear();
    recv_size.clear();
    recv_from.clear();

    std::size_t TotalRcvsVolume = 0;
    for (const auto& kv : m_RcvVols) // loop over senders
    {
        std::size_t nbytes = 0;
        for (const auto& cct : kv.second)
        {
            nbytes += (*this)[cct.dstIndex].nBytes(cct.dbox,icomp,ncomp);
        }

        BL_ASSERT(nbytes < std::numeric_limits<int>::max());

        TotalRcvsVolume += nbytes;

        recv_data.push_back(nullptr);
        recv_size.push_back(static_cast<int>(nbytes));
        recv_from.push_back(kv.first);
    }

    if (TotalRcvsVolume == 0)
    {
        the_recv_data = nullptr;
    }
    else
    {
        the_recv_data = static_cast<char*>(BLPgas::alloc(TotalRcvsVolume));

        int nrecv = recv_from.size();

        std::size_t offset = 0;
        for (int i = 0; i < nrecv; ++i)
        {
            if (recv_size[i] > 0)
            {
                recv_data[i] = the_recv_data + offset;

                /// Send an AM to the sender with the recv pointer
                upcxx::global_ptr<void> dst_ptr =
                    upcxx::global_ptr<void>(recv_data[i], upcxx::myrank());
                // Increment the event reference before launching the remote task
                recv_event->incref();
                // Launch a remote task on the sender to send me data
                BLPgas::Request(recv_from[i], dst_ptr, recv_size[i], SeqNum, recv_event);
                upcxx::advance(); // poll the UPC++ progress engine and the network
                offset += recv_size[i];
            }
        }
    }
}
#endif /* BL_USE_UPCXX */


template <class FAB>
template <class>
void
FabArray<FAB>::BuildMask (const Box& phys_domain, const Periodicity& period,
			  value_type covered, value_type notcovered, 
			  value_type physbnd, value_type interior)
{
    int ncomp = this->nComp();
    int ngrow = this->nGrow();

    const FabArrayBase::FB& TheFB = this->getFB(period);

    const CopyComTagsContainer&      LocTags = *(TheFB.m_LocTags);
    const MapOfCopyComTagContainers& RcvTags = *(TheFB.m_RcvTags);

    Box domain = amrex::convert(phys_domain, boxArray().ixType());
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	if (period.isPeriodic(i)) {
	    domain.grow(i, ngrow);
	}
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this, true); mfi.isValid(); ++mfi)
    {
	FAB& fab = (*this)[mfi];

	Box gbx = mfi.growntilebox();
	fab.setVal(physbnd, gbx, 0, ncomp);

	gbx &= domain;
	fab.setVal(notcovered, gbx, 0, ncomp);

	const Box& tbx = mfi.tilebox();
	fab.setVal(interior, tbx, 0, ncomp);
    }

    int N_locs = LocTags.size();
#ifdef _OPENMP
#pragma omp parallel for if (TheFB.m_threadsafe_loc)
#endif
    for (int i = 0; i < N_locs; ++i) {
	const CopyComTag& tag = LocTags[i];
	(*this)[tag.dstIndex].setVal(covered, tag.dbox, 0, ncomp);
    }

    for (MapOfCopyComTagContainers::const_iterator it = RcvTags.begin(); it != RcvTags.end(); ++it) {
	int N = it->second.size();
#ifdef _OPENMP
#pragma omp parallel for if (TheFB.m_threadsafe_rcv)
#endif
	for (int i = 0; i < N; ++i) {
	    const CopyComTag& tag = it->second[i];
	    (*this)[tag.dstIndex].setVal(covered, tag.dbox, 0, ncomp);
	}
    }
}

template <typename FAB> std::map<int, std::map<int, FabArray<FAB> *> > FabArray<FAB>::allocatedFAPointers;
template <typename FAB> std::set<int> FabArray<FAB>::noallocFAPIds;

}

#endif /*BL_FABARRAY_H*/
