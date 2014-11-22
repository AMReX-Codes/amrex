
#ifndef BL_FABARRAY_H
#define BL_FABARRAY_H

#include <iostream>
#include <cstring>
#include <limits>
#include <map>
#include <utility>
#include <vector>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <BLassert.H>
#include <Array.H>

#include <Box.H>
#include <BoxLib.H>
#include <BoxArray.H>
#include <BoxDomain.H> 
#include <FArrayBox.H>
#include <DistributionMapping.H>
#include <ParallelDescriptor.H>
#include <ccse-mpi.H>
#include <Profiler.H>

//
// Helper class
//

class FillBoxId
{
  public:

    FillBoxId ()
        :
        m_fillBoxId(-1),
        m_fabIndex(-1)
        {}
    FillBoxId (int newid, const Box& fillbox)
        :
        m_fillBox(fillbox),
        m_fillBoxId(newid),
        m_fabIndex(-1)
        {}

    int Id () const              { return m_fillBoxId;    }
    int FabIndex () const        { return m_fabIndex;     }
    void FabIndex (int fabindex) { m_fabIndex = fabindex; }
    const Box& box () const      { return m_fillBox;      }

private:

    Box m_fillBox;
    int m_fillBoxId;
    int m_fabIndex;
};

//
// This is meant to be a concrete class not a polymorphic one.
//

class FabArrayBase
{
public:

    FabArrayBase ();

    ~FabArrayBase();
    //
    // Returns the grow factor that defines the region of definition.
    //
    int nGrow () const { return n_grow; }
    //
    // Returns number of variables associated with each point (nvar).
    //
    int nComp () const { return n_comp; }
    //
    // Returns a constant reference to the BoxArray that defines the
    // valid region associated with this FabArray.
    //
    const BoxArray& boxArray () const { return boxarray; }
    //
    // Returns a constant reference to the Kth Box in the BoxArray.
    // That is, the valid region of the Kth grid.
    //
    const Box& box (int K) const { return boxarray[K]; }
    //
    // Returns the Kth FABs Box in the FabArray.
    // That is, the region the Kth fab is actually defined on.
    //
    const Box fabbox (int K) const;
    //
    // Returns the number of FABs in the FabArray..
    //
    int size () const { return boxarray.size(); }
    //
    // Returns constant reference to associated DistributionMapping.
    //
    const DistributionMapping& DistributionMap () const { return distributionMap; }
    //
    // Returns constant reference to indices in the FabArray that we own.
    //
    const Array<int>& IndexMap () const { return indexMap; }
    //
    // Returns local index in the vector of FABs.
    //
    int localindex (int K) const { 
	std::vector<int>::const_iterator it = std::find(indexMap.begin(),indexMap.end(),K);
	if (it != indexMap.end()) {
	    return it - indexMap.begin();
	}
	else {
	    return -1;
	}
    }
    //
    // Flush the cache of self-intersection info used by FillBoundary.
    //
    static void FlushSICache ();
    //
    // The size of the cache of self-intersection info.
    //
    static int SICacheSize ();
    //
    // Some static member templates used throughout the code.
    //
    template<typename T>
    static void SetRecvTag (std::map< int,std::vector<T> >& m_RcvTags,
                            int                             src_owner,
                            const T&                        tag,
                            std::map<int,int>&              m_RcvVols,
                            const Box&                      bx);

    template<typename T>
    static void SetSendTag (std::map< int,std::vector<T> >& m_SndTags,
                            int                             dst_owner,
                            const T&                        tag,
                            std::map<int,int>&              m_SndVols,
                            const Box&                      bx);

    template<typename T>
    static void GrokAsyncSends (int                 N_snds,
                                Array<MPI_Request>& send_reqs,
                                Array<T*>&          send_data,
                                Array<MPI_Status>&  stats);

    template<typename T1, typename T2>
    static void PostRcvs (const std::map< int,std::vector<T1> >& m_RcvTags,
                          const std::map<int,int>&               m_RcvVols,
                          T2*&                                   the_recv_data,
                          Array<T2*>&                            recv_data,
                          Array<int>&                            recv_from,
                          Array<MPI_Request>&                    recv_reqs,
                          int                                    ncomp,
                          int                                    SeqNum);
    //
    // Used by a bunch of routines when communicating via MPI.
    //
    struct CopyComTag
    {
        Box box;
        int fabIndex;
        int srcIndex;
        //
        // Some typedefs & helper functions used throughout the code.
        //
        typedef std::vector<CopyComTag> CopyComTagsContainer;

        typedef std::map<int,CopyComTagsContainer> MapOfCopyComTagContainers;
    };
    //
    // Used in caching self-intersection info for FillBoundary().
    //
    struct SI
    {
        SI ();

        SI (const BoxArray&            ba,
            const DistributionMapping& dm,
            int                        ngrow,
            bool                       cross);

        ~SI ();

        bool operator== (const SI& rhs) const;
        bool operator!= (const SI& rhs) const { return !operator==(rhs); }

        int bytes () const;
        //
        // Some useful typedefs.
        //
        typedef CopyComTag::CopyComTagsContainer CopyComTagsContainer;

        typedef CopyComTag::MapOfCopyComTagContainers MapOfCopyComTagContainers;
        //
        // Basic data.
        //
        BoxArray            m_ba;
        DistributionMapping m_dm;
        int                 m_ngrow;
        bool                m_cross;
        bool                m_reused;
        //
        // The cache of local and send/recv per FillBoundary().
        //
        CopyComTagsContainer*      m_LocTags;
        MapOfCopyComTagContainers* m_SndTags;
        MapOfCopyComTagContainers* m_RcvTags;
        std::map<int,int>*         m_SndVols;
        std::map<int,int>*         m_RcvVols;
    };
    //
    // Some useful typedefs for the FillBoundary() cache.
    //
    typedef std::multimap<int,FabArrayBase::SI> FBCache;

    typedef FBCache::iterator FBCacheIter;

    static FBCache m_TheFBCache;
    //
    // When copy()ing from one FabArray to another we can do a copy or an add.
    //
    enum CpOp { COPY = 0, ADD = 1 };

    struct CPC
    {
        CPC ();

        CPC (const BoxArray&            dstba,
             const BoxArray&            srcba,
             const DistributionMapping& dstdm,
             const DistributionMapping& srcdm);

        ~CPC ();

        bool operator== (const CPC& rhs) const;
        bool operator!= (const CPC& rhs) const { return !operator==(rhs); }

        int bytes () const;

        static void FlushCache ();

        BoxArray            m_dstba;
        BoxArray            m_srcba;
        DistributionMapping m_dstdm;
        DistributionMapping m_srcdm;
        bool                m_reused;
        //
        // Some useful typedefs.
        //
        typedef CopyComTag::CopyComTagsContainer CopyComTagsContainer;

        typedef CopyComTag::MapOfCopyComTagContainers MapOfCopyComTagContainers;
        //
        // The cache of local and send/recv info per FabArray::copy().
        //
        CopyComTagsContainer*      m_LocTags;
        MapOfCopyComTagContainers* m_SndTags;
        MapOfCopyComTagContainers* m_RcvTags;
        std::map<int,int>*         m_SndVols;
        std::map<int,int>*         m_RcvVols;
    };
    //
    // Some useful typedefs for the copy() cache.
    //
    typedef std::multimap<int,FabArrayBase::CPC> CPCCache;

    typedef CPCCache::iterator CPCCacheIter;

    static CPCCache m_TheCopyCache;

    static CPCCacheIter TheCPC (const CPC&          cpc,
                                const FabArrayBase& dst,
                                const FabArrayBase& src);

    static void EraseFromTheCPC (CPCCacheIter& cache_it);
    //
    // Used for collecting information used in communicating FABs.
    //
    struct FabComTag
    {
        int fromProc;
        int toProc;
        int fabIndex;
        int fineIndex;
        int srcComp;
        int destComp;
        int nComp;
        int face;
        int fabArrayId;
        int fillBoxId;
        int procThatNeedsData;
        int procThatHasData;
        Box box;

        FabComTag ()
            :
            fromProc(0),
            toProc(0),
            fabIndex(0),
            fineIndex(0),
            srcComp(0),
            destComp(0),
            nComp(0),
            face(0),
            fabArrayId(0),
            fillBoxId(0),
            procThatNeedsData(0),
            procThatHasData(0) {}
    };
    //
    // Returns cached self-intersection records or builds them.
    //
    static FBCacheIter TheFB (bool cross, const FabArrayBase& mf);
    //
    // Default tilesize in MFIter
    //
    static IntVect mfiter_tile_size;
    //
    // The maximum number of components to copy() at a time.
    //
    static int MaxComp;
    //
    // Use MPI_Asend() instead of MPI_Send() in CollectData() and copy().
    //
    // Turn on via ParmParse using "fabarray.do_async_sends=1" in inputs file.
    //
    // Default is false.
    //
    static bool do_async_sends;
    //
    // Shuffle the order in which we send MPI data in copy(), FillBoundary() and
    // FillPeriodicBoundary() in hopes of evening out any networking bottlenecks.
    //
    // Default is false.
    //
    static bool do_random_shuffle;
    //
    // Print out some stuff; default is false.
    //
    static bool Verbose;
    //
    // Initialize from ParmParse with "fabarray" prefix.
    //
    static void Initialize ();
    static void Finalize ();

protected:
    //
    // The data ...
    //
    mutable BoxArray    boxarray;   /* So FabSets can modify'm */
    DistributionMapping distributionMap;
    Array<int>          indexMap;
    int                 n_grow;
    int                 n_comp;
};

class MFIter
{
public:
    //
    // Construct a MFIter.
    //
    explicit MFIter (const FabArrayBase& fabarray)
        :
        fabArray(fabarray),
        currentIndex(0),
	tileSize(D_DECL(10240,10240,10240))
	{
	    Initialize();
	}
    explicit MFIter (const FabArrayBase& fabarray, bool do_tiling)
        :
        fabArray(fabarray),
        currentIndex(0),
	tileSize(D_DECL(10240,10240,10240))
	{
	    if (do_tiling) 
		tileSize = FabArrayBase::mfiter_tile_size;
	    Initialize();
	}
    explicit MFIter (const FabArrayBase& fabarray, const IntVect& tilesize)
        :
        fabArray(fabarray),
        currentIndex(0),
	tileSize(tilesize)
	{
	    Initialize();
	}
    //
    // Returns the tile Box at the current index.
    //
    const Box& tilebox () const { return tileArray[currentIndex]; }
    //
    // Returns the valid Box that current tile resides.
    //
    const Box& validbox () const { return fabArray.box(indexMap[currentIndex]); }
    //
    // Returns the Box of the FAB at which we currently point.
    //
    const Box fabbox () const { return fabArray.fabbox(indexMap[currentIndex]); }
    //
    // Increments iterator to the next tile we own.
    //
    void operator++ () { ++currentIndex;}
    //
    // Is the iterator valid i.e. is it associated with a FAB?
    //
    bool isValid () { return currentIndex < indexMap.size(); }
    //
    // The index into the underlying BoxArray of the current FAB.
    //
    int index () const { return indexMap[currentIndex]; }
    //
    // local index
    //
    int CurrentIndex () const { return currentIndex; }
    //
    // Constant reference to FabArray over which we're iterating.
    //
    const FabArrayBase& theFabArrayBase () const { return fabArray; }

private:

    const FabArrayBase& fabArray;
    int                 currentIndex;
    IntVect             tileSize;
    Array<int>          indexMap;
    Array<Box>          tileArray;

    void Initialize () {

	int tid = 0;
	int nthreads = 1;

#ifdef _OPENMP
	if (omp_in_parallel()) {
	    tid = omp_get_thread_num();
	    nthreads = omp_get_num_threads();
	}
#endif

	// First let's figure out how many tiles we are going to have
	int n_tot_tiles=0;
	Array<IntVect> nt_in_fab;
	for (int i=0; i<fabArray.IndexMap().size(); i++) {
	    int K = fabArray.IndexMap()[i]; 
	    const Box bx = fabArray.box(K);

	    int ntiles = 1;
	    IntVect nt;
	    for (int d=0; d<BL_SPACEDIM; d++) {
		nt[d] = std::max(bx.length(d)/tileSize[d], 1);
		ntiles *= nt[d];
	    }

	    nt_in_fab.push_back(nt);
	    n_tot_tiles += ntiles;
	}

	// the tile no range for this thread
	int tlo, thi;
	if (n_tot_tiles < nthreads) { // there are more threads than tiles
	    if (tid < n_tot_tiles) {
		tlo = thi = tid;
	    }
	    else {
		return;
	    }
	}
	else {
	    int tiles_per_thread = n_tot_tiles/nthreads;
	    int nleft = n_tot_tiles - tiles_per_thread*nthreads;
	    if (tid < nleft) {
		tlo = tid*(tiles_per_thread+1);
		thi = tlo + tiles_per_thread;
	    }
	    else {
		tlo = tid*tiles_per_thread + nleft;
		thi = tlo + tiles_per_thread - 1;
	    }
	}

	indexMap.reserve(thi-tlo+1);
	tileArray.reserve(thi-tlo+1);

	int it = 0;
	for (int i=0; i<fabArray.IndexMap().size(); i++) {

	    int ntiles = 1;
	    for (int d=0; d<BL_SPACEDIM; d++) {
		ntiles *= nt_in_fab[i][d];
	    }

	    if (it+ntiles-1 < tlo) {
		it += ntiles;
		continue;
	    }

	    int K = fabArray.IndexMap()[i]; 
	    const Box bx = fabArray.box(K);

	    IntVect tsize, nleft;
	    for (int d=0; d<BL_SPACEDIM; d++) {
		int ncells = bx.length(d);
		tsize[d] = ncells/nt_in_fab[i][d];
		nleft[d] = ncells - nt_in_fab[i][d]*tsize[d];
	    }

	    IntVect small, big, ijk;  // note that the initial values are all zero.
	    ijk[0] = -1;
	    for (int t=0; t<ntiles; t++) {
		for (int d=0; d<BL_SPACEDIM; d++) {
		    if (ijk[d]<nt_in_fab[i][d]-1) {
			ijk[d]++;
			break;
		    }
		    else {
			ijk[d] = 0;
		    }
		}
		
		for (int d=0; d<BL_SPACEDIM; d++) {
		    if (ijk[d] < nleft[d]) {
			small[d] = ijk[d]*(tsize[d]+1);
			big[d] = small[d] + tsize[d];
		    }
		    else {
			small[d] = ijk[d]*tsize[d] + nleft[d];
			big[d] = small[d] + tsize[d] - 1;
		    }
		}

		indexMap.push_back(K);
		Box tbx(small, big, bx.ixType());
		tbx.shift(bx.smallEnd());
		tileArray.push_back(tbx);
		it++;
		if (it > thi) return;
	    }
	}
    }
};

//
// A forward declaration.
//
template <class FAB> class FabArray;
template <class FAB> class FabArrayCopyDescriptor;

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
// An enumumeration that controls whether or not the memory for a FAB
// will actually be allocated on construction of a FabArray.
// Possible values are: Fab_noallocate and Fab_allocate.
//

enum FabAlloc { Fab_noallocate = 0, Fab_allocate };

template <class FAB>
class FabArray
    :
    public FabArrayBase
{
public:

    typedef typename FAB::value_type value_type;
    //
    // Constructs an empty FabArray<FAB>.
    //
    FabArray ();
    //
    // Construct a FabArray<FAB> with a valid region defined by bxs
    // and a region of definition defined by the grow factor ngrow
    // and the number of components nvar.
    // If mem_mode is defined to be Fab_allocate then FABs are
    // allocated for each Box in the BoxArray.  The size of the Kth
    // FAB is given by bxs[K] grown by ngrow.  If mem_mode is defined
    // to be Fab_noallocate, then no FABs are allocated at this time,
    // but can be defined later.  The number of components in each
    // FAB is not specified and is expected to be implicit in the
    // definition of the FAB class.  That is, the FAB constructor will
    // take only a Box argument.  Call this constructor number two.
    //
    FabArray (const BoxArray& bxs,
              int             nvar,
              int             ngrow,
              FabAlloc        mem_mode = Fab_allocate);

    FabArray (const BoxArray&            bxs,
              int                        nvar,
              int                        ngrow,
              const DistributionMapping& dm,
              FabAlloc                   mem_mode = Fab_allocate);
    //
    // The destructor -- deletes all FABs in the array.
    //
    ~FabArray ();
    //
    // Define this FabArray identically to that performed by
    // the constructor having an analogous function signature.
    // This is only valid if this FabArray was defined using
    // the default constructor.
    //
    void define (const BoxArray& bxs,
                 int             nvar,
                 int             ngrow,
                 FabAlloc        mem_mode);

    void define (const BoxArray&            bxs,
                 int                        nvar,
                 int                        ngrow,
		 const DistributionMapping& dm,
                 FabAlloc                   mem_mode);
    //
    // Returns true if the FabArray is well-defined.  That is,
    // if FABs are allocated for each Box in the BoxArray and the
    // sizes of the FABs and the number of components are consistent
    // with the definition of the FabArray.
    //
    bool ok () const;
    //
    // Returns a constant reference to the FAB associated with the Kth element.
    //
    const FAB& operator[] (const MFIter& mfi) const;

    const FAB& get (const MFIter& mfi) const { return this->operator[](mfi); }
    //
    // Returns a reference to the FAB associated mfi.
    //
    FAB& operator[] (const MFIter& mfi);

    FAB& get (const MFIter& mfi) { return this->operator[](mfi); }
    //
    // Returns a constant reference to the FAB associated with the Kth element.
    //
    const FAB& operator[] (int K) const;

    const FAB& get (int K) const { return this->operator[](K); }
    //
    // Returns a reference to the FAB associated with the Kth element.
    //
    FAB& operator[] (int K);

    FAB& get (int K)  { return this->operator[](K); }
    //
    // Explicitly set the Kth FAB in the FabArray to point to elem.
    //
    void setFab (int K, FAB* elem);

    void setFab (const MFIter&mfi, FAB* elem);
    //
    // Releases FAB memory in the FabArray.
    //
    void clear ();
    //
    // Set all components in the entire region of each FAB to val.
    //
    void setVal (value_type val);
    void operator= (const value_type& val);
    //
    // Set the value of num_comp components in the valid region of
    // each FAB in the FabArray, starting at component comp to val.
    // Also set the value of nghost boundary cells.
    //
    void setVal (value_type val,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);
    //
    // Set the value of num_comp components in the valid region of
    // each FAB in the FabArray, starting at component comp, as well
    // as nghost boundary cells, to val, provided they also intersect
    // with the Box region.
    //
    void setVal (value_type val,
                 const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);
    //
    // Set all components in the valid region of each FAB in the
    // FabArray to val, including nghost boundary cells.
    //
    void setVal (value_type val,
                 int        nghost);
    //
    // Set all components in the valid region of each FAB in the
    // FabArray to val, including nghost boundary cells, that also
    // intersect the Box region.
    //
    void setVal (value_type val,
                 const Box& region,
                 int        nghost);
    //
    // Set all values in the boundary region to val.
    //
    void setBndry (value_type val);
    //
    // Set ncomp values in the boundary region, starting at start_comp to val.
    //
    void setBndry (value_type val,
                   int        strt_comp,
                   int        ncomp);
    //
    // This function copies data from fa to this FabArray.  Each FAB
    // in fa is intersected with all FABs in this FabArray and a copy
    // is performed on the region of intersection.  The intersection
    // is restricted to the valid region of each FAB.
    //
    void copy (const FabArray<FAB>& fa,
               CpOp                 op = FabArrayBase::COPY);
    //
    // This function copies data from src to this FabArray.  Each FAB
    // in src is intersected with all FABs in this FabArray and a copy
    // is performed on the region of intersection.  The intersection
    // is restricted to the num_comp components starting at src_comp
    // in the FabArray src, with the destination components in this
    // FabArray starting at dest_comp.
    //
    void copy (const FabArray<FAB>& src,
               int                  src_comp,
               int                  dest_comp,
               int                  num_comp,
               CpOp                 op = FabArrayBase::COPY);
    //
    // Copies the values contained in the intersection of the
    // valid region of this FabArray with the FAB dest into dest.
    //
    void copy (FAB& dest) const;
    //
    // Copies the values contained in the intersection of the
    // valid region of this FabArray with the FAB dest and the Box
    // subbox into that subregion of dest.
    //
    void copy (FAB&       dest,
               const Box& subbox) const;
    //
    // Copies the values contained in the intersection of the
    // num_comp component valid region of this FabArray, starting at
    // component src_comp, with the FAB dest into dest, starting at
    // component dest_comp in dest.
    //
    void copy (FAB& dest,
               int  src_comp,
               int  dest_comp,
               int  num_comp) const;
    //
    // Copies the values contained in the intersection of the
    // num_comp component valid region of this FabArray, starting at
    // component src_comp, with the FAB dest and the Box subbox, into
    // dest, starting at component dest_comp in dest.
    //
    void copy (FAB&       dest,
               const Box& subbox,
               int        src_comp,
               int        dest_comp,
               int        num_comp) const;

    void shift (const IntVect& v);

    bool defined (int i) const;
    bool defined (const MFIter& mfi) const;
    //
    // Copy on intersection within a FabArray.  Data is copied from
    // valid regions to intersecting regions of definition.  The
    // purpose is to fill in the boundary regions of each FAB in
    // the FabArray.
    //
    void FillBoundary (bool cross = false);
    //
    // Same as FillBoundary(), but only copies ncomp components starting at scomp.
    //
    void FillBoundary (int scomp, int ncomp, bool cross = false);

protected:
    //
    // Helper function for define().
    //
    void defineDoit (const BoxArray&            bxs,
                     int                        nvar,
                     int                        ngrow,
                     FabAlloc                   mem_mode,
                     const DistributionMapping* dm);
    //
    // The data.
    //
    std::vector<FAB*> m_fabs;

private:
    typedef typename std::vector<FAB*>::iterator    Iterator;
    //
    // These are disallowed.
    //
    FabArray (const FabArray<FAB>&);
    FabArray<FAB>& operator= (const FabArray<FAB>&);
    //
    // This is used locally in all define functions.
    //
    void AllocFabs ();
};

class FabArrayId
{
public:

    explicit FabArrayId (int newid = -1)
        :
        fabArrayId(newid) {}

    int Id () const { return fabArrayId; }

    bool operator== (const FabArrayId& rhs) const
    {
        return fabArrayId == rhs.fabArrayId;
    }

private:

    int fabArrayId;
};

//
// This enum and the FabCopyDescriptor class should really be nested
// in FabArrayCopyDescriptor (not done for portability reasons).
//

enum FillType { FillLocally, FillRemotely, Unfillable };

template <class FAB>
struct FabCopyDescriptor
{
    FabCopyDescriptor ();

    ~FabCopyDescriptor ();

    FAB*     localFabSource;
    Box      subBox;
    int      myProc;
    int      copyFromProc;
    int      copyFromIndex;
    int      fillBoxId;
    int      srcComp;
    int      destComp;
    int      nComp;
    FillType fillType;
    bool     cacheDataAllocated;

private:
    //
    // Disallowed.
    //
    FabCopyDescriptor (const FabCopyDescriptor&);
    FabCopyDescriptor& operator= (const FabCopyDescriptor&);
};

template <class FAB>
FabCopyDescriptor<FAB>::FabCopyDescriptor ()
    :
    localFabSource(0),
    myProc(-1),
    copyFromProc(-1),
    copyFromIndex(-1),
    fillBoxId(-1),
    srcComp(-1),
    destComp(-1),
    nComp(-1),
    fillType(Unfillable),
    cacheDataAllocated(false)
{}

template <class FAB>
FabCopyDescriptor<FAB>::~FabCopyDescriptor ()
{
    if (cacheDataAllocated)
        delete localFabSource;
}

//
// This class orchestrates filling a destination fab of size destFabBox
// from fabarray on the local processor (myProc).
//

template <class FAB>
class FabArrayCopyDescriptor
{
  typedef std::multimap<int,FabCopyDescriptor<FAB>*> FCDMap;

  typedef typename FCDMap::value_type     FCDMapValueType;
  typedef typename FCDMap::iterator       FCDMapIter;
  typedef typename FCDMap::const_iterator FCDMapConstIter;

  public:

    FabArrayCopyDescriptor ();

    ~FabArrayCopyDescriptor ();

    FabArrayId RegisterFabArray(FabArray<FAB> *fabarray);

    FillBoxId AddBox (FabArrayId fabarrayid,
                      const Box& destFabBox,
                      BoxList*   unfilledBoxes);

    FillBoxId AddBox (FabArrayId fabarrayid,
                      const Box& destFabBox,
                      BoxList*   unfilledBoxes,
                      int        srccomp,
                      int        destcomp,
                      int        numcomp);
    //
    // Add a box but only from FabArray[fabarrayindex].
    //
    FillBoxId AddBox (FabArrayId fabarrayid,
                      const Box& destFabBox,
                      BoxList*   unfilledBoxes,
                      int        fabarrayindex,
                      int        srccomp,
                      int        destcomp,
                      int        numcomp,
                      bool       bUseValidBox = true);

    void CollectData ();

    void FillFab (FabArrayId       fabarrayid,
                  const FillBoxId& fillboxid,
                  FAB&             destFab);

    void FillFab (FabArrayId       fabarrayid,
                  const FillBoxId& fillboxid,
                  FAB&             destFab,
                  const Box&       destBox);

    void PrintStats () const;

    bool DataAvailable () const { return dataAvailable; }

    void clear ();

    int nFabArrays () const { return fabArrays.size(); }

    int nFabComTags () const { return fabComTagList.size(); }

    int nFabCopyDescs () const { return fabCopyDescList.size(); }

private:
    //
    // These are disallowed.
    //
    FabArrayCopyDescriptor (const FabArrayCopyDescriptor<FAB>&);

    FabArrayCopyDescriptor<FAB>& operator= (const FabArrayCopyDescriptor<FAB> &);
    //
    // Helper function for AddBox() routines.
    //
    void AddBoxDoIt (FabArrayId fabarrayid,
                     const Box& destFabBox,
                     BoxList*   returnedUnfilledBoxes,
                     int        faindex,
                     int        srccomp,
                     int        destcomp,
                     int        numcomp,
                     bool       bUseValidBox,
                     BoxDomain& unfilledBoxDomain);
    //
    // Some useful typedefs.
    //
    typedef std::map<int,int> IntIntMap;

    typedef std::vector<FabArrayBase::FabComTag> FabComTagContainer;

    typedef std::vector<FabComTagContainer::const_iterator> FabComTagIterContainer;
    //
    // The data.
    //
    std::vector<FabArray<FAB>*> fabArrays;
    std::vector<FCDMap>         fabCopyDescList;
    FabComTagContainer          fabComTagList;
    int                         nextFillBoxId;
    bool                        dataAvailable;
};

template<typename T>
void
FabArrayBase::SetRecvTag (std::map< int,std::vector<T> >& m_RcvTags,
                          int                             src_owner,
                          const T&                        tag,
                          std::map<int,int>&              m_RcvVols,
                          const Box&                      bx)
{
    m_RcvTags[src_owner].push_back(tag);

    std::map<int,int>::iterator vol_it = m_RcvVols.find(src_owner);

    const int vol = bx.numPts();

    if (vol_it != m_RcvVols.end())
    {
        vol_it->second += vol;
    }
    else
    {
        m_RcvVols[src_owner] = vol;
    }
}

template<typename T>
void
FabArrayBase::SetSendTag (std::map< int,std::vector<T> >& m_SndTags,
                          int                             dst_owner,
                          const T&                        tag,
                          std::map<int,int>&              m_SndVols,
                          const Box&                      bx)
{
    m_SndTags[dst_owner].push_back(tag);

    std::map<int,int>::iterator vol_it = m_SndVols.find(dst_owner);

    const int vol = bx.numPts();

    if (vol_it != m_SndVols.end())
    {
        vol_it->second += vol;
    }
    else
    {
        m_SndVols[dst_owner] = vol;
    }
}

template<typename T1, typename T2>
void
FabArrayBase::PostRcvs (const std::map< int,std::vector<T1> >& m_RcvTags,
                        const std::map<int,int>&               m_RcvVols,
                        T2*&                                   the_recv_data,
                        Array<T2*>&                            recv_data,
                        Array<int>&                            recv_from,
                        Array<MPI_Request>&                    recv_reqs,
                        int                                    ncomp,
                        int                                    SeqNum)
{
    int TotalRcvsVolume = 0;

    for (std::map<int,int>::const_iterator it = m_RcvVols.begin(),
             End = m_RcvVols.end();
         it != End;
         ++it)
    {
        TotalRcvsVolume += it->second;
    }

    TotalRcvsVolume *= ncomp;

    BL_ASSERT((TotalRcvsVolume*sizeof(T2)) < std::numeric_limits<int>::max());

    the_recv_data = static_cast<T2*>(BoxLib::The_Arena()->alloc(TotalRcvsVolume*sizeof(T2)));

    int Offset = 0;

    for (typename std::map< int,std::vector<T1> >::const_iterator m_it = m_RcvTags.begin(),
             m_End = m_RcvTags.end();
         m_it != m_End;
         ++m_it)
    {
        std::map<int,int>::const_iterator vol_it = m_RcvVols.find(m_it->first);

        BL_ASSERT(vol_it != m_RcvVols.end());

        const int N = vol_it->second*ncomp;

        BL_ASSERT(N < std::numeric_limits<int>::max());

        recv_data.push_back(&the_recv_data[Offset]);
        recv_from.push_back(m_it->first);
        recv_reqs.push_back(ParallelDescriptor::Arecv(recv_data.back(),N,m_it->first,SeqNum).req());

        Offset += N;
    }
}

template<typename T>
void
FabArrayBase::GrokAsyncSends (int                 N_snds,
                              Array<MPI_Request>& send_reqs,
                              Array<T*>&          send_data,
                              Array<MPI_Status>&  stats)
{
#ifdef BL_USE_MPI
    BL_ASSERT(FabArrayBase::do_async_sends && N_snds > 0);

    stats.resize(N_snds);

    BL_ASSERT(send_reqs.size() == N_snds);
    BL_ASSERT(send_data.size() == N_snds);

    Array<int> indx;
    BL_COMM_PROFILE_WAITSOME(Profiler::Waitall, send_reqs, N_snds, indx, stats, false);

    BL_MPI_REQUIRE( MPI_Waitall(N_snds, send_reqs.dataPtr(), stats.dataPtr()) );

    BL_COMM_PROFILE_WAITSOME(Profiler::Waitall, send_reqs, N_snds, indx, stats, false);

    for (int i = 0; i < N_snds; i++)
        BoxLib::The_Arena()->free(send_data[i]);
#endif /*BL_USE_MPI*/
}

template <class FAB>
bool
FabArray<FAB>::defined (int K) const
{
    int li = localindex(K);
    if (li >= 0 && li < m_fabs.size() && m_fabs[li] != 0) {
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
    int li = mfi.CurrentIndex();
    if (li < m_fabs.size() && m_fabs[li] != 0) {
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
    BL_ASSERT(mfi.CurrentIndex() < indexMap.size());
    return *m_fabs[mfi.CurrentIndex()];
}

template <class FAB>
FAB&
FabArray<FAB>::operator[] (const MFIter& mfi)
{
    BL_ASSERT(mfi.CurrentIndex() < indexMap.size());
    return *m_fabs[mfi.CurrentIndex()];
}

template <class FAB>
const FAB&
FabArray<FAB>::operator[] (int K) const
{
    int li = localindex(K);
    BL_ASSERT(li >=0 && li < indexMap.size());
    return *m_fabs[li];
}

template <class FAB>
FAB&
FabArray<FAB>::operator[] (int K)
{
    int li = localindex(K);
    BL_ASSERT(li >=0 && li < indexMap.size());
    return *m_fabs[li];
}

template <class FAB>
void
FabArray<FAB>::clear ()
{
    for (Iterator it = m_fabs.begin(); it != m_fabs.end(); ++it) 
	delete *it;
    
    m_fabs.clear();
    boxarray.clear();
}

template <class FAB>
void
FabArray<FAB>::setVal (value_type val,
                       int        nghost)
{
    setVal(val,0,n_comp,nghost);
}

template <class FAB>
void
FabArray<FAB>::setVal (value_type   val,
                         const Box& region,
                         int        nghost)
{
    setVal(val,region,0,n_comp,nghost);
}

template <class FAB>
FabArray<FAB>::FabArray ()
{}

template <class FAB>
FabArray<FAB>::FabArray (const BoxArray& bxs,
                         int             nvar,
                         int             ngrow,
                         FabAlloc        alloc)
{
    define(bxs,nvar,ngrow,alloc);
}

template <class FAB>
FabArray<FAB>::FabArray (const BoxArray&            bxs,
                         int                        nvar,
                         int                        ngrow,
                         const DistributionMapping& dm,
                         FabAlloc                   alloc)
{
    define(bxs,nvar,ngrow,dm,alloc);
}

template <class FAB>
FabArray<FAB>::~FabArray ()
{
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
            if (get(fai).box() != BoxLib::grow(box(fai.index()),n_grow))
            {
                isok = false;
            }
        }
        else
        {
            isok = false;
        }
    }

    ParallelDescriptor::ReduceLongAnd(isok);

    return isok != 0;
}

template <class FAB>
void
FabArray<FAB>::defineDoit (const BoxArray&            bxs,
                           int                        nvar,
                           int                        ngrow,
                           FabAlloc                   alloc,
                           const DistributionMapping* dm)
{
    BL_ASSERT(ngrow >= 0);
    BL_ASSERT(boxarray.size() == 0);

    n_grow = ngrow;
    n_comp = nvar;

    boxarray.define(bxs);

    if (dm == 0)
    {
        distributionMap.define(boxarray,ParallelDescriptor::NProcs());
    }
    else
    {
        BL_ASSERT(dm->ProcessorMap().size() == bxs.size()+1);

        distributionMap = *dm;
    }

    const int MyProc = ParallelDescriptor::MyProc();

    for (int i = 0, N = boxarray.size(); i < N; i++)
        if (distributionMap[i] == MyProc)
            indexMap.push_back(i);

    if (alloc == Fab_allocate)
        AllocFabs();
}

template <class FAB>
void
FabArray<FAB>::define (const BoxArray& bxs,
                       int             nvar,
                       int             ngrow,
                       FabAlloc        alloc)
{
    defineDoit(bxs,nvar,ngrow,alloc,0);
}

template <class FAB>
void
FabArray<FAB>::define (const BoxArray&            bxs,
                       int                        nvar,
                       int                        ngrow,
                       const DistributionMapping& dm,
                       FabAlloc                   alloc)
{
    defineDoit(bxs,nvar,ngrow,alloc,&dm);
}

template <class FAB>
void
FabArray<FAB>::AllocFabs ()
{
    m_fabs.reserve(indexMap.size());

    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        const Box tmp = BoxLib::grow(fai.validbox(), n_grow);
	m_fabs.push_back(new FAB(tmp, n_comp));
    }
}

template <class FAB>
void
FabArray<FAB>::setFab (int  boxno,
                       FAB* elem)
{
    //
    // Must check it is of the proper size.
    //
    if (n_comp == 0)
        n_comp = elem->nComp();

    BL_ASSERT(n_comp == elem->nComp());
    BL_ASSERT(boxarray.size() > 0);
    BL_ASSERT(elem->box() == BoxLib::grow(boxarray[boxno],n_grow));
    BL_ASSERT(!this->defined(boxno));
    BL_ASSERT(distributionMap[boxno] == ParallelDescriptor::MyProc());

    if (m_fabs.size() == 0) {
	m_fabs.resize(indexMap.size());
    }

    m_fabs[localindex(boxno)] = elem;
}

template <class FAB>
void
FabArray<FAB>::setFab (const MFIter& mfi,
                       FAB* elem)
{
    //
    // Must check it is of the proper size.
    //
    if (n_comp == 0)
        n_comp = elem->nComp();

    BL_ASSERT(n_comp == elem->nComp());
    BL_ASSERT(boxarray.size() > 0);
    BL_ASSERT(elem->box() == BoxLib::grow(boxarray[mfi.index()],n_grow));
    BL_ASSERT(!this->defined(mfi));
    BL_ASSERT(distributionMap[mfi.index()] == ParallelDescriptor::MyProc());

    if (m_fabs.size() == 0) {
	m_fabs.resize(indexMap.size());
    }

    m_fabs[mfi.CurrentIndex()] = elem;
}

template <class FAB>
void
FabArray<FAB>::setBndry (value_type val)
{
    setBndry(val, 0, n_comp);
}

template <class FAB>
void
FabArray<FAB>::setBndry (value_type val,
                         int        strt_comp,
                         int        ncomp)
{
    if (n_grow > 0)
    {
        for (MFIter fai(*this); fai.isValid(); ++fai)
        {
            get(fai).setComplement(val, fai.validbox(), strt_comp, ncomp);
        }
    }
}

template <class FAB>
void
FabArray<FAB>::copy (const FabArray<FAB>& src,
                     int                  scomp,
                     int                  dcomp,
                     int                  ncomp,
                     CpOp                 op)
{
    BL_PROFILE("FabArray::copy()");

    if (size() == 0 || src.size() == 0) return;

    BL_ASSERT(op == FabArrayBase::COPY || op == FabArrayBase::ADD);
    BL_ASSERT(boxArray()[0].ixType() == src.boxArray()[0].ixType());

    if ((src.boxArray()[0].cellCentered() || op == FabArrayBase::COPY) &&
        (boxarray == src.boxarray && distributionMap == src.distributionMap))
    {
        //
        // Short-circuit full intersection code if we're doing copy()s or if
        // we're doing plus()s on cell-centered data.  Don't do plus()s on
        // non-cell-centered data this simplistic way.
        //
        for (MFIter fai(*this); fai.isValid(); ++fai)
        {
            const Box& bx = fai.validbox();

            if (op == FabArrayBase::COPY)
            {
                get(fai).copy(src[fai],bx,scomp,bx,dcomp,ncomp);
            }
            else
            {
                get(fai).plus(src[fai],bx,bx,scomp,dcomp,ncomp);
            }
        }

        return;
    }

    const CPC cpc(boxarray, src.boxarray, distributionMap, src.distributionMap);

    FabArrayBase::CPCCacheIter cache_it = FabArrayBase::TheCPC(cpc, *this, src);

    BL_ASSERT(cache_it != FabArrayBase::m_TheCopyCache.end());

    const CPC& thecpc = cache_it->second;

    if (ParallelDescriptor::NProcs() == 1)
    {
        //
        // There can only be local work to do.
        //
        for (CPC::CopyComTagsContainer::const_iterator it = thecpc.m_LocTags->begin(),
                 End = thecpc.m_LocTags->end();
             it != End;
             ++it)
        {
            const CopyComTag& tag = *it;

            if (op == FabArrayBase::COPY)
            {
                get(tag.fabIndex).copy(src[tag.srcIndex],tag.box,scomp,tag.box,dcomp,ncomp);
            }
            else
            {
                get(tag.fabIndex).plus(src[tag.srcIndex],tag.box,tag.box,scomp,dcomp,ncomp);
            }
        }

        return;
    }

#ifdef BL_USE_MPI
    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    const int SeqNum = ParallelDescriptor::SeqNum();

    if (thecpc.m_LocTags->empty() && thecpc.m_RcvTags->empty() && thecpc.m_SndTags->empty())
        //
        // No work to do.
        //
        return;
    //
    // Send/Recv at most MaxComp components at a time to cut down memory usage.
    //
    int NCompLeft = ncomp;

    for (int i = 0, SC = scomp, DC = dcomp; i < ncomp; )
    {
        const int NC = std::min(NCompLeft,FabArrayBase::MaxComp);

        Array<MPI_Status>  stats;
        Array<int>         recv_from;
        Array<value_type*> recv_data, send_data;
        Array<MPI_Request> recv_reqs, send_reqs;
        //
        // Post rcvs. Allocate one chunk of space to hold'm all.
        //
        value_type* the_recv_data = 0;

        FabArrayBase::PostRcvs(*thecpc.m_RcvTags,*thecpc.m_RcvVols,the_recv_data,recv_data,recv_from,recv_reqs,NC,SeqNum);
        //
        // Send the FAB data.
        //
        FAB fab;
        //
        // Get the indices of SndTags and shuffle'm a bit.
        //
        Array<int> indices;

        for (CPC::MapOfCopyComTagContainers::const_iterator m_it = thecpc.m_SndTags->begin(),
                 m_End = thecpc.m_SndTags->end();
             m_it != m_End;
             ++m_it)
        {
            indices.push_back(m_it->first);
        }

        if (FabArrayBase::do_random_shuffle)
            std::random_shuffle(indices.begin(), indices.end());

        for (int j = 0, NN = indices.size(); j < NN; j++)
        {
            CPC::MapOfCopyComTagContainers::const_iterator m_it = thecpc.m_SndTags->find(indices[j]);

            std::map<int,int>::const_iterator vol_it = thecpc.m_SndVols->find(m_it->first);

            BL_ASSERT(vol_it != thecpc.m_SndVols->end());

            const int N = vol_it->second*NC;

            BL_ASSERT(N < std::numeric_limits<int>::max());

            value_type* data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(N*sizeof(value_type)));
            value_type* dptr = data;

            for (CPC::CopyComTagsContainer::const_iterator it = m_it->second.begin(),
                     End = m_it->second.end();
                 it != End;
                 ++it)
            {
                const Box& bx = it->box;
                src[it->srcIndex].copyToMem(bx,SC,NC,dptr);
                const int Cnt = bx.numPts()*NC;
                dptr += Cnt;
            }
            BL_ASSERT(data+N == dptr);

            if (FabArrayBase::do_async_sends)
            {
                send_data.push_back(data);
                send_reqs.push_back(ParallelDescriptor::Asend(data,N,m_it->first,SeqNum).req());
            }
            else
            {
                ParallelDescriptor::Send(data,N,m_it->first,SeqNum);
                BoxLib::The_Arena()->free(data);
            }
        }
        //
        // Do the local work.  Hope for a bit of communication/computation overlap.
        //
        for (CPC::CopyComTagsContainer::const_iterator it = thecpc.m_LocTags->begin(),
                 End = thecpc.m_LocTags->end();
             it != End;
             ++it)
        {
            const CopyComTag& tag = *it;

            BL_ASSERT(    distributionMap[tag.fabIndex] == ParallelDescriptor::MyProc());
            BL_ASSERT(src.distributionMap[tag.srcIndex] == ParallelDescriptor::MyProc());

            if (op == FabArrayBase::COPY)
            {
                get(tag.fabIndex).copy(src[tag.srcIndex],tag.box,SC,tag.box,DC,NC);
            }
            else
            {
                get(tag.fabIndex).plus(src[tag.srcIndex],tag.box,tag.box,SC,DC,NC);
            }
        }

	//
	//  wait and unpack
	//

        const int N_rcvs = thecpc.m_RcvTags->size();

	if (N_rcvs > 0)
	{
	    stats.resize(N_rcvs);
	    
	    BL_MPI_REQUIRE( MPI_Waitall(N_rcvs, recv_reqs.dataPtr(), stats.dataPtr()) );
	    
	    for (int k = 0; k < N_rcvs; k++)
	    {
		const value_type* dptr = recv_data[k];
		
		BL_ASSERT(dptr != 0);
		
		CPC::MapOfCopyComTagContainers::const_iterator m_it = thecpc.m_RcvTags->find(recv_from[k]);
		
		BL_ASSERT(m_it != thecpc.m_RcvTags->end());
		
		for (CPC::CopyComTagsContainer::const_iterator it = m_it->second.begin(),
			 End = m_it->second.end();
		     it != End;
		     ++it)
		{
		    const Box& bx  = it->box;
		    const int  Cnt = bx.numPts()*NC;
		    
		    if (op == FabArrayBase::COPY)
		    {
			get(it->fabIndex).copyFromMem(bx,DC,NC,dptr);
		    }
		    else
		    {
			fab.resize(bx,NC);
			memcpy(fab.dataPtr(), dptr, Cnt*sizeof(value_type));
			
			get(it->fabIndex).plus(fab,bx,bx,0,DC,NC);
		    }
		    
		    dptr += Cnt;
		}
	    }
	}
	
        BoxLib::The_Arena()->free(the_recv_data);

        if (FabArrayBase::do_async_sends && !thecpc.m_SndTags->empty())
            FabArrayBase::GrokAsyncSends(thecpc.m_SndTags->size(),send_reqs,send_data,stats);

        i         += NC;
        SC        += NC;
        DC        += NC;
        NCompLeft -= NC;
    }

#endif /*BL_USE_MPI*/
}

template <class FAB>
void
FabArray<FAB>::copy (const FabArray<FAB>& src, CpOp op)
{
    copy(src,0,0,nComp(),op);
}

//
// Copies to FABs, note that destination is first arg.
//

template <class FAB>
void
FabArray<FAB>::copy (FAB& dest) const
{
    copy(dest, dest.box(), 0, 0, dest.nComp());
}

template <class FAB>
void
FabArray<FAB>::copy (FAB&       dest,
                     const Box& subbox) const
{
    copy(dest, subbox, 0, 0, dest.nComp());
}

template <class FAB>
void
FabArray<FAB>::copy (FAB& dest,
                     int  scomp,
                     int  dcomp,
                     int  ncomp) const
{
    copy(dest, dest.box(), scomp, dcomp, ncomp);
}

template <class FAB>
void
FabArray<FAB>::copy (FAB&       dest,
                     const Box& subbox,
                     int        scomp,
                     int        dcomp,
                     int        ncomp) const
{
    BL_ASSERT(dcomp + ncomp <= dest.nComp());

    if (ParallelDescriptor::NProcs() == 1)
    {
        for (int j = 0, N = size(); j < N; ++j)
        {
            if (boxarray[j].intersects(subbox))
            {
                Box destbox = boxarray[j] & subbox;

                dest.copy(get(j),destbox,scomp,destbox,dcomp,ncomp);
            }
        }

        return;
    }

    FAB ovlp;

    for (int i = 0, N = size(); i < N; i++)
    {
        if (subbox.intersects(boxarray[i]))
        {
            Box bx = subbox & boxarray[i];

            ovlp.resize(bx,ncomp);

            if (ParallelDescriptor::MyProc() == distributionMap[i])
            {
                ovlp.copy(get(i),bx,scomp,bx,0,ncomp);
            }

            const int N = bx.numPts()*ncomp;

            ParallelDescriptor::Bcast(ovlp.dataPtr(),N,distributionMap[i]);

            dest.copy(ovlp,bx,0,bx,dcomp,ncomp);
        }
    }
}

template <class FAB>
void
FabArray<FAB>::setVal (value_type val)
{
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        get(fai).setVal(val);
    }
}

template <class FAB>
void
FabArray<FAB>::operator= (const value_type& val)
{
    setVal(val);
}

template <class FAB>
void
FabArray<FAB>::setVal (value_type val,
                       int        comp,
                       int        ncomp,
                       int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+ncomp <= n_comp);

    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        get(fai).setVal(val,BoxLib::grow(fai.validbox(),nghost), comp, ncomp);
    }
}

template <class FAB>
void
FabArray<FAB>::setVal (value_type val,
                       const Box& region,
                       int        comp,
                       int        ncomp,
                       int        nghost)
{
    BL_ASSERT(nghost >= 0 && nghost <= n_grow);
    BL_ASSERT(comp+ncomp <= n_comp);

    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        Box b = BoxLib::grow(fai.validbox(),nghost) & region;

        if (b.ok())
            get(fai).setVal(val, b, comp, ncomp);
    }
}

template <class FAB>
void
FabArray<FAB>::shift (const IntVect& v)
{
    for(int id(0); id < BL_SPACEDIM; ++id)
    {
      boxarray.shift(id, v[id]);
    }
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        get(fai).shift(v);
    }
}

template <class FAB>
FabArrayCopyDescriptor<FAB>::FabArrayCopyDescriptor ()
    :
    nextFillBoxId(0),
    dataAvailable(false)
{}

template <class FAB>
FabArrayId
FabArrayCopyDescriptor<FAB>::RegisterFabArray(FabArray<FAB>* fabarray)
{
    BL_ASSERT(fabArrays.size() == fabCopyDescList.size());

    FabArrayId result(fabArrays.size());

    fabArrays.push_back(fabarray);  /* Bump size() by one */

    fabCopyDescList.push_back(FCDMap());

    return result;
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::AddBoxDoIt (FabArrayId fabarrayid,
                                         const Box& destFabBox,
                                         BoxList*   returnedUnfilledBoxes,
                                         int        faindex,
                                         int        srccomp,
                                         int        destcomp,
                                         int        numcomp,
                                         bool       bUseValidBox,
                                         BoxDomain& unfilledBoxDomain)
{
    const int MyProc = ParallelDescriptor::MyProc();

    FabArray<FAB>* fabArray = fabArrays[fabarrayid.Id()];

    BL_ASSERT(faindex >= 0 && faindex < fabArray->size());

    Box intersect = destFabBox;

    if (bUseValidBox)
    {
        intersect &= fabArray->box(faindex);
    }
    else
    {
        intersect &= fabArray->fabbox(faindex);
    }

    if (intersect.ok())
    {
        FabCopyDescriptor<FAB>* fcd = new FabCopyDescriptor<FAB>;

        int remoteProc     = fabArray->DistributionMap()[faindex];
        fcd->fillBoxId     = nextFillBoxId;
        fcd->subBox        = intersect;
        fcd->myProc        = MyProc;
        fcd->copyFromProc  = remoteProc;
        fcd->copyFromIndex = faindex;
        fcd->srcComp       = srccomp;
        fcd->destComp      = destcomp;
        fcd->nComp         = numcomp;

        if (MyProc == remoteProc)
        {
            //
            // Data is local.
            //
            fcd->fillType       = FillLocally;
            fcd->localFabSource = &(*fabArray)[faindex];
        }
        else
        {
            //
            // Data is remote.
            //
            FabArrayBase::FabComTag fabComTag;

            dataAvailable               = false;
            fcd->fillType               = FillRemotely;
            fcd->localFabSource         = new FAB(intersect, numcomp);
            fcd->cacheDataAllocated     = true;
            fabComTag.fabArrayId        = fabarrayid.Id();
            fabComTag.fillBoxId         = nextFillBoxId;
            fabComTag.fabIndex          = faindex;
            fabComTag.procThatNeedsData = MyProc;
            fabComTag.procThatHasData   = remoteProc;
            fabComTag.box               = intersect;
            fabComTag.srcComp           = srccomp;
            fabComTag.destComp          = destcomp;
            fabComTag.nComp             = numcomp;
            //
            // Do not send the data yet.
            //
            fabComTagList.push_back(fabComTag);
        }

        fabCopyDescList[fabarrayid.Id()].insert(FCDMapValueType(fcd->fillBoxId,fcd));

        if (returnedUnfilledBoxes != 0)
        {
            unfilledBoxDomain.rmBox(intersect);
        }
    }
}

template <class FAB>
FillBoxId
FabArrayCopyDescriptor<FAB>::AddBox (FabArrayId fabarrayid,
                                     const Box& destFabBox,
                                     BoxList*   returnedUnfilledBoxes,
                                     int        srccomp,
                                     int        destcomp,
                                     int        numcomp)
{
    BoxDomain unfilledBoxDomain(destFabBox.ixType());

    if (returnedUnfilledBoxes != 0)
    {
        unfilledBoxDomain.add(destFabBox);
    }

    std::vector< std::pair<int,Box> > isects;

    fabArrays[fabarrayid.Id()]->boxArray().intersections(destFabBox,isects);

    for (int j = 0, N = isects.size(); j < N; j++)
    {
        AddBoxDoIt(fabarrayid,
                   destFabBox,
                   returnedUnfilledBoxes,
                   isects[j].first,
                   srccomp,
                   destcomp,
                   numcomp,
                   true,
                   unfilledBoxDomain);
    }

    if (returnedUnfilledBoxes != 0)
    {
        returnedUnfilledBoxes->clear();
        (*returnedUnfilledBoxes) = unfilledBoxDomain.boxList();
    }

    return FillBoxId(nextFillBoxId++, destFabBox);
}

template <class FAB>
FillBoxId
FabArrayCopyDescriptor<FAB>::AddBox (FabArrayId fabarrayid,
                                     const Box& destFabBox,
                                     BoxList*   returnedUnfilledBoxes,
                                     int        fabarrayindex,
                                     int        srccomp,
                                     int        destcomp,
                                     int        numcomp,
                                     bool       bUseValidBox)
{
    BoxDomain unfilledBoxDomain(destFabBox.ixType());

    if (returnedUnfilledBoxes != 0)
    {
        unfilledBoxDomain.add(destFabBox);
    }

    AddBoxDoIt(fabarrayid,
               destFabBox,
               returnedUnfilledBoxes,
               fabarrayindex,
               srccomp,
               destcomp,
               numcomp,
               bUseValidBox,
               unfilledBoxDomain);

    if (returnedUnfilledBoxes != 0)
    {
        returnedUnfilledBoxes->clear();
        (*returnedUnfilledBoxes) = unfilledBoxDomain.boxList();
    }

    return FillBoxId(nextFillBoxId++, destFabBox);
}

template <class FAB>
FillBoxId
FabArrayCopyDescriptor<FAB>::AddBox (FabArrayId fabarrayid,
                                     const Box& destFabBox,
                                     BoxList*   returnedUnfilledBoxes)
{
    return AddBox(fabarrayid,
                  destFabBox,
                  returnedUnfilledBoxes,
                  0,
                  0,
                  fabArrays[fabarrayid.Id()]->nComp(),
                  true);
}

template <class FAB>
FabArrayCopyDescriptor<FAB>::~FabArrayCopyDescriptor()
{
   clear();
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::clear ()
{
    for (unsigned int i = 0, N = fabCopyDescList.size(); i < N; ++i)
    {
        for (FCDMapIter fmi = fabCopyDescList[i].begin(), End = fabCopyDescList[i].end();
             fmi != End;
             ++fmi)
        {
            delete (*fmi).second;
        }
    }

    fabArrays.clear();
    fabCopyDescList.clear();
    fabComTagList.clear();

    nextFillBoxId = 0;
    dataAvailable = false;
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::CollectData ()
{
    dataAvailable = true;

    if (ParallelDescriptor::NProcs() == 1) return;

#if BL_USE_MPI
    typedef typename FAB::value_type value_type;
    //
    // Make sure we can treat CommData as a stream of integers.
    //
    BL_ASSERT(sizeof(ParallelDescriptor::CommData) == ParallelDescriptor::CommData::DIM*sizeof(int));

    BL_PROFILE("FabArrayCopyDescriptor::CollectData()");

    const int MyProc = ParallelDescriptor::MyProc();

    int Total_Rcvs_Size = 0;
    //
    // We use this to make finding matching FabComTags more efficient.
    //
    std::map< int,FabComTagIterContainer > RcvTags;

    IntIntMap Snds, Rcvs, Npts;
    //
    // Set Rcvs[i] to # of blocks needed from CPU i
    //
    for (FabComTagContainer::const_iterator it = fabComTagList.begin(),
             End = fabComTagList.end();
         it != End;
         ++it)
    {
        BL_ASSERT(it->box.ok());
        BL_ASSERT(it->procThatNeedsData == MyProc);
        BL_ASSERT(it->procThatHasData   != MyProc);

        const int Who = it->procThatHasData;
        const int Cnt = (it->box.numPts())*(it->nComp);

        RcvTags[Who].push_back(it);

        Total_Rcvs_Size += Cnt;

        if (Rcvs.count(Who) > 0)
        {
            Rcvs[Who] += 1;
        }
        else
        {
            Rcvs[Who] = 1;
        }

        if (Npts.count(Who) > 0)
        {
            Npts[Who] += Cnt;
        }
        else
        {
            Npts[Who] = Cnt;
        }
    }
    BL_ASSERT(Rcvs.count(MyProc) == 0);

    BL_ASSERT((Total_Rcvs_Size*sizeof(value_type)) < std::numeric_limits<int>::max());

    const int NProcs = ParallelDescriptor::NProcs();

    {
        Array<int> SndsArray(NProcs,0), RcvsArray(NProcs,0);

        for (IntIntMap::const_iterator it = Rcvs.begin(), End = Rcvs.end(); it != End; ++it)
            RcvsArray[it->first] = it->second;

        {
            BL_PROFILE("CollectData_Alltoall()");
	    BL_COMM_PROFILE(Profiler::Alltoall, sizeof(int), ParallelDescriptor::MyProc(),
	                    Profiler::BeforeCall());

            BL_MPI_REQUIRE( MPI_Alltoall(RcvsArray.dataPtr(),
                                         1,
                                         ParallelDescriptor::Mpi_typemap<int>::type(),
                                         SndsArray.dataPtr(),
                                         1,
                                         ParallelDescriptor::Mpi_typemap<int>::type(),
                                         ParallelDescriptor::Communicator()) );

	    BL_COMM_PROFILE(Profiler::Alltoall, sizeof(int), ParallelDescriptor::MyProc(),
	                    Profiler::AfterCall());

        }
        BL_ASSERT(SndsArray[MyProc] == 0);

        for (int i = 0; i < NProcs; i++)
            if (SndsArray[i] > 0)
                Snds[i] = SndsArray[i];
    }

    Array<ParallelDescriptor::CommData> cd_others_need, cd_that_i_need;

    {
        Array<int> sdispls(NProcs,0), rdispls(NProcs,0), scnts(NProcs,0), rcnts(NProcs,0);

        int nrcvs = 0;
        for (IntIntMap::const_iterator it = Snds.begin(), End = Snds.end(); it != End; ++it)
        {
            nrcvs           += it->second;
            rcnts[it->first] = it->second;
        }
        for (int i = 1; i < NProcs; i++)
            rdispls[i] = rdispls[i-1] + rcnts[i-1];

        int nsnds = 0;
        for (IntIntMap::const_iterator it = Rcvs.begin(), End = Rcvs.end(); it != End; ++it)
        {
            nsnds           += it->second;
            scnts[it->first] = it->second;
        }
        for (int i = 1; i < NProcs; i++)
            sdispls[i] = sdispls[i-1] + scnts[i-1];

        Array<int> index(sdispls);

        cd_others_need.resize(nrcvs+1); // +1 so there's always at least one element.
        cd_that_i_need.resize(nsnds+1); // +1 so there's always at least one element.

        for (FabComTagContainer::const_iterator it = fabComTagList.begin(),
                 End = fabComTagList.end();
             it != End;
             ++it)
        {
            ParallelDescriptor::CommData data(0,
                                              it->fabIndex,
                                              MyProc,
                                              0,
                                              it->nComp,
                                              it->srcComp,
                                              it->fabArrayId,
                                              it->box);

            cd_that_i_need[index[it->procThatHasData]++] = data;
        }
        //
        // Increment displacements to indicate integers not CommData.
        //
        for (int i = 0; i < NProcs; i++)   scnts[i] *= ParallelDescriptor::CommData::DIM;
        for (int i = 0; i < NProcs; i++)   rcnts[i] *= ParallelDescriptor::CommData::DIM;
        for (int i = 1; i < NProcs; i++) sdispls[i] *= ParallelDescriptor::CommData::DIM;
        for (int i = 1; i < NProcs; i++) rdispls[i] *= ParallelDescriptor::CommData::DIM;

        {
            BL_PROFILE("CollectData_Alltoallv()");
	    BL_COMM_PROFILE(Profiler::Alltoallv, nrcvs * sizeof(int), ParallelDescriptor::MyProc(),
	                    Profiler::BeforeCall());

            BL_MPI_REQUIRE( MPI_Alltoallv(cd_that_i_need.dataPtr(),
                                          scnts.dataPtr(),
                                          sdispls.dataPtr(),
                                          ParallelDescriptor::Mpi_typemap<int>::type(),
                                          cd_others_need.dataPtr(),
                                          rcnts.dataPtr(),
                                          rdispls.dataPtr(),
                                          ParallelDescriptor::Mpi_typemap<int>::type(),
                                          ParallelDescriptor::Communicator()) );

	    BL_COMM_PROFILE(Profiler::Alltoallv, nrcvs * sizeof(int), ParallelDescriptor::MyProc(),
	                    Profiler::AfterCall());
        }
        cd_that_i_need.clear();
    }

    Array<int>         roffset, who_R;
    Array<MPI_Status>  stats;
    Array<MPI_Request> recv_reqs, send_reqs;
    Array<value_type*> send_data;
    //
    // Post receives.  Allocate data for rcvs as one big chunk.
    //
    const int   SeqNum    = ParallelDescriptor::SeqNum();
    value_type* recv_data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(Total_Rcvs_Size*sizeof(value_type)));

    int Idx = 0;
    for (IntIntMap::const_iterator it = Rcvs.begin(), End = Rcvs.end(); it != End; ++it)
    {
        const int Who = it->first;
        const int Cnt = Npts[Who];
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        who_R.push_back(Who);
        recv_reqs.push_back(ParallelDescriptor::Arecv(&recv_data[Idx],Cnt,Who,SeqNum).req());
        roffset.push_back(Idx);
        Idx += Cnt;
    }
    //
    // Send the FAB data.
    //
    FAB fab;

    Idx = 0;
    for (IntIntMap::const_iterator it = Snds.begin(), End = Snds.end(); it != End; ++it)
    {
        const int Who   = it->first;
        const int NSnds = it->second;

        int N = 0;
        const ParallelDescriptor::CommData* cd = &cd_others_need[Idx];
        for (int k = 0; k < NSnds; k++, cd++)
            N += cd->box().numPts()*cd->nComp();

        BL_ASSERT(N < std::numeric_limits<int>::max());

        value_type* data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(N*sizeof(value_type)));
        value_type* dptr = data;

        cd = &cd_others_need[Idx];

        for (int k = 0; k < NSnds; k++, cd++)
        {
            BL_ASSERT(cd->id() == 0);
            BL_ASSERT(cd->fromproc() == Who);
            const int Cnt = cd->box().numPts()*cd->nComp();
            (*fabArrays[cd->fabarrayid()])[cd->fabindex()].copyToMem(cd->box(),cd->srcComp(),cd->nComp(),dptr);
            dptr += Cnt;
        }

        BL_ASSERT(data+N == dptr);

        if (FabArrayBase::do_async_sends)
        {
            send_data.push_back(data);
            send_reqs.push_back(ParallelDescriptor::Asend(data,N,Who,SeqNum).req());
        }
        else
        {
            ParallelDescriptor::Send(data,N,Who,SeqNum);
            BoxLib::The_Arena()->free(data);
        }

        Idx += NSnds;
    }

    cd_others_need.clear();

    //
    //  wait and unpack
    //

    int N_rcvs = recv_reqs.size();

    if (N_rcvs > 0) 
    {
	stats.resize(N_rcvs);
	
	std::pair<FCDMapIter,FCDMapIter> match;
	std::map< int,FabComTagIterContainer >::const_iterator found;
	
	BL_MPI_REQUIRE( MPI_Waitall(N_rcvs, recv_reqs.dataPtr(), stats.dataPtr()) );
	
	for (int k = 0; k < N_rcvs; k++)
	{
	    const int         Who     = who_R[k];
	    const value_type* dptr    = &recv_data[roffset[k]];
	    
	    BL_ASSERT(dptr != 0);
	    
	    found = RcvTags.find(Who);
	    
	    BL_ASSERT(found != RcvTags.end());
	    
	    const FabComTagIterContainer& tags = found->second;
	    
	    for (FabComTagIterContainer::const_iterator it = tags.begin(), End = tags.end();
		 it != End;
		 ++it)
	    {
		const FabArrayBase::FabComTag& tag = **it;                  
		
		BL_ASSERT(tag.procThatHasData == Who);
		
		match = fabCopyDescList[tag.fabArrayId].equal_range(tag.fillBoxId);
		
		for (FCDMapIter fmi = match.first; fmi != match.second; ++fmi)
		{
		    FabCopyDescriptor<FAB>* fcdp = (*fmi).second;
		    
		    BL_ASSERT(fcdp->fillBoxId == tag.fillBoxId);
		    
		    if (fcdp->subBox == tag.box)
		    {
			BL_ASSERT(fcdp->localFabSource->dataPtr() != 0);
			BL_ASSERT(fcdp->localFabSource->box() == tag.box);
			const int Cnt = tag.box.numPts()*tag.nComp;
			fcdp->localFabSource->copyFromMem(tag.box,0,tag.nComp,dptr);
			dptr += Cnt;
			break;
		    }
		}
	    }
	}
    }

    BoxLib::The_Arena()->free(recv_data);

    if (FabArrayBase::do_async_sends && !send_reqs.empty())
    {
        //
        // Now grok the asynchronous send buffers & free up send buffer space.
        //
        const int N_snds = send_reqs.size();

        stats.resize(N_snds);

        BL_COMM_PROFILE(Profiler::Waitall, sizeof(value_type), Profiler::BeforeCall(), N_snds);
        BL_MPI_REQUIRE( MPI_Waitall(N_snds, send_reqs.dataPtr(), stats.dataPtr()) );
        BL_COMM_PROFILE(Profiler::Waitall, sizeof(value_type), Profiler::AfterCall(), N_snds);

        for (int i = 0; i < N_snds; i++)
            BoxLib::The_Arena()->free(send_data[i]);
    }

#endif /*BL_USE_MPI*/
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::FillFab (FabArrayId       faid,
                                      const FillBoxId& fillboxid,
                                      FAB&             destFab)
{
    BL_ASSERT(dataAvailable);

    std::pair<FCDMapIter,FCDMapIter> match = fabCopyDescList[faid.Id()].equal_range(fillboxid.Id());

    for (FCDMapIter fmi = match.first; fmi != match.second; ++fmi)
    {
        FabCopyDescriptor<FAB>* fcdp = (*fmi).second;

        BL_ASSERT(fcdp->fillBoxId == fillboxid.Id());

        destFab.copy(*fcdp->localFabSource,
                     fcdp->subBox,
                     fcdp->fillType == FillLocally ? fcdp->srcComp : 0,
                     fcdp->subBox,
                     fcdp->destComp,
                     fcdp->nComp);
    }
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::FillFab (FabArrayId       faid,
                                        const FillBoxId& fillboxid,
                                        FAB&             destFab,
                                        const Box&       destBox)
{
    BL_ASSERT(dataAvailable);

    FCDMapIter fmi = fabCopyDescList[faid.Id()].lower_bound(fillboxid.Id());

    BL_ASSERT(fmi != fabCopyDescList[faid.Id()].end());

    FabCopyDescriptor<FAB>* fcdp = (*fmi).second;

    BL_ASSERT(fcdp->fillBoxId == fillboxid.Id());

    BL_ASSERT(fcdp->subBox.sameSize(destBox));

    destFab.copy(*fcdp->localFabSource,
                 fcdp->subBox,
                 fcdp->fillType == FillLocally ? fcdp->srcComp : 0,
                 destBox,
                 fcdp->destComp,
                 fcdp->nComp);

    BL_ASSERT(++fmi == fabCopyDescList[faid.Id()].upper_bound(fillboxid.Id()));
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::PrintStats () const
{
    const int MyProc = ParallelDescriptor::MyProc();

    std::cout << "----- "
         << MyProc
         << ":  Parallel stats for FabArrayCopyDescriptor:" << '\n';

    for (int fa = 0; fa < fabArrays.size(); ++fa)
    {
      std::cout << "fabArrays["
             << fa
             << "]->boxArray() = "
             << fabArrays[fa]->boxArray()
             << '\n';
    }
}

template <class FAB>
void
FabArray<FAB>::FillBoundary (bool cross)
{
    FillBoundary(0, nComp(), cross);
}

template <class FAB>
void
FabArray<FAB>::FillBoundary (int  scomp,
                             int  ncomp,
                             bool cross)
{
    BL_PROFILE("FabArray::FillBoundary()");

    if ( n_grow <= 0 ) return;

    FabArrayBase::FBCacheIter cache_it = FabArrayBase::TheFB(cross,*this);

    BL_ASSERT(cache_it != FabArrayBase::m_TheFBCache.end());

    const FabArrayBase::SI& TheSI = cache_it->second;

    if (ParallelDescriptor::NProcs() == 1)
    {
        //
        // There can only be local work to do.
        //
        for (SI::CopyComTagsContainer::const_iterator it = TheSI.m_LocTags->begin(),
                 End = TheSI.m_LocTags->end();
             it != End;
             ++it)
        {
            const CopyComTag& tag = *it;

            BL_ASSERT(distributionMap[tag.fabIndex] == ParallelDescriptor::MyProc());
            BL_ASSERT(distributionMap[tag.srcIndex] == ParallelDescriptor::MyProc());

            get(tag.fabIndex).copy(get(tag.srcIndex),tag.box,scomp,tag.box,scomp,ncomp);
        }

        return;
    }

#ifdef BL_USE_MPI
    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    const int SeqNum = ParallelDescriptor::SeqNum();

    if (TheSI.m_LocTags->empty() && TheSI.m_RcvTags->empty() && TheSI.m_SndTags->empty())
        //
        // No work to do.
        //
        return;

    Array<MPI_Status>  stats;
    Array<int>         recv_from;
    Array<value_type*> recv_data, send_data;
    Array<MPI_Request> recv_reqs, send_reqs;
    //
    // Post rcvs. Allocate one chunk of space to hold'm all.
    //
    value_type* the_recv_data = 0;

    FabArrayBase::PostRcvs(*TheSI.m_RcvTags,*TheSI.m_RcvVols,the_recv_data,recv_data,recv_from,recv_reqs,ncomp,SeqNum);
    //
    // Send the FAB data.
    //
    FAB fab;
    //
    // Get the indices of SndTags and shuffle'm a bit.
    //
    Array<int> indices;

    for (SI::MapOfCopyComTagContainers::const_iterator m_it = TheSI.m_SndTags->begin(),
             m_End = TheSI.m_SndTags->end();
         m_it != m_End;
         ++m_it)
    {
        indices.push_back(m_it->first);
    }

    if (FabArrayBase::do_random_shuffle)
        std::random_shuffle(indices.begin(), indices.end());

    for (int j = 0, NN = indices.size(); j < NN; j++)
    {
        SI::MapOfCopyComTagContainers::const_iterator m_it = TheSI.m_SndTags->find(indices[j]);

        std::map<int,int>::const_iterator vol_it = TheSI.m_SndVols->find(m_it->first);

        BL_ASSERT(vol_it != TheSI.m_SndVols->end());

        const int N = vol_it->second*ncomp;

        BL_ASSERT(N < std::numeric_limits<int>::max());

        value_type* data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(N*sizeof(value_type)));
        value_type* dptr = data;

        for (SI::CopyComTagsContainer::const_iterator it = m_it->second.begin(),
                 End = m_it->second.end();
             it != End;
             ++it)
        {
            BL_ASSERT(distributionMap[it->srcIndex] == ParallelDescriptor::MyProc());
            const Box& bx = it->box;
            get(it->srcIndex). copyToMem(bx,scomp,ncomp,dptr);
            const int Cnt = bx.numPts()*ncomp;
            dptr += Cnt;
        }
        BL_ASSERT(data+N == dptr);

        if (FabArrayBase::do_async_sends)
        {
            send_data.push_back(data);
            send_reqs.push_back(ParallelDescriptor::Asend(data,N,m_it->first,SeqNum).req());
        }
        else
        {
            ParallelDescriptor::Send(data,N,m_it->first,SeqNum);
            BoxLib::The_Arena()->free(data);
        }
    }
    //
    // Do the local work.  Hope for a bit of communication/computation overlap.
    //
    for (SI::CopyComTagsContainer::const_iterator it = TheSI.m_LocTags->begin(),
             End = TheSI.m_LocTags->end();
         it != End;
         ++it)
    {
        const CopyComTag& tag = *it;

        BL_ASSERT(distributionMap[tag.fabIndex] == ParallelDescriptor::MyProc());
        BL_ASSERT(distributionMap[tag.srcIndex] == ParallelDescriptor::MyProc());

        get(tag.fabIndex).copy(get(tag.srcIndex),tag.box,scomp,tag.box,scomp,ncomp);
    }

    //
    //  wait and unpack
    //

    const int N_rcvs = TheSI.m_RcvTags->size();

    if (N_rcvs > 0)
    {
	stats.resize(N_rcvs);
	
	BL_MPI_REQUIRE( MPI_Waitall(N_rcvs, recv_reqs.dataPtr(), stats.dataPtr()) );
	
	for (int k = 0; k < N_rcvs; k++)
	{
	    const value_type* dptr = recv_data[k];
	    
	    BL_ASSERT(dptr != 0);
	    
	    SI::MapOfCopyComTagContainers::const_iterator m_it = TheSI.m_RcvTags->find(recv_from[k]);
	    
	    BL_ASSERT(m_it != TheSI.m_RcvTags->end());
	    
	    for (SI::CopyComTagsContainer::const_iterator it = m_it->second.begin(),
		     End = m_it->second.end();
		 it != End;
		 ++it)
	    {
		const Box& bx  = it->box;
		const int  Cnt = bx.numPts()*ncomp;
		get(it->fabIndex).copyFromMem(bx,scomp,ncomp,dptr);
		dptr += Cnt;
	    }
	}
    }

    BoxLib::The_Arena()->free(the_recv_data);

    if (FabArrayBase::do_async_sends && !TheSI.m_SndTags->empty())
        FabArrayBase::GrokAsyncSends(TheSI.m_SndTags->size(),send_reqs,send_data,stats);

#endif /*BL_USE_MPI*/
}

#endif /*BL_FABARRAY_H*/
