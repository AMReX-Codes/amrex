
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
#include <BLProfiler.H>

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
    bool empty () const { return boxarray.empty(); }
    //
    // Returns a constant reference to the BoxArray that defines the
    // valid region associated with this FabArray.
    //
    const BoxArray& boxArray () const { return boxarray; }
    //
    // Returns the Kth Box in the BoxArray.
    // That is, the valid region of the Kth grid.
    //
    Box box (int K) const { return boxarray[K]; }
    //
    // Returns the Kth FABs Box in the FabArray.
    // That is, the region the Kth fab is actually defined on.
    //
    Box fabbox (int K) const;
    //
    // Returns the number of FABs in the FabArray..
    //
    int size () const { return boxarray.size(); }
    //
    // Returns the number of local FABs in the FabArray..
    //
    int local_size () const { return indexMap.size(); }
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
	std::vector<int>::const_iterator low
	    = std::lower_bound(indexMap.begin(), indexMap.end(), K);
	if (low != indexMap.end() && *low == K) {
	    return low - indexMap.begin();
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

    template<typename T>
    static void PostRcvs (const std::map<int,int>&               m_RcvVols,
                          T*&                                    the_recv_data,
                          Array<T*>&                             recv_data,
                          Array<int>&                            recv_from,
                          Array<MPI_Request>&                    recv_reqs,
                          int                                    ncomp,
                          int                                    SeqNum);
    //
    struct CacheStats
    {
	int         size;     // current size: nbuild - nerase
	int         maxsize;  // highest water mark of size
	int         maxuse;   // max # of uses of a cached item
	long        nuse;     // # of uses of the whole cache
	long        nbuild;   // # of build operations
	long        nerase;   // # of erase operations
	std::string name;     // name of the cache
	CacheStats (const std::string& name_) 
	    : size(0),maxsize(0),maxuse(0),nuse(0),nbuild(0),nerase(0),name(name_) {;}
	void recordBuild () {
	    ++size;  
	    ++nbuild;  
	    maxsize = std::max(maxsize, size); 
	}
	void recordErase (int n) { 
	    // n: how many times the item to be deleted has been used.
	    --size;
	    ++nerase;
	    maxuse = std::max(maxuse, n);
	}
	void recordUse () { ++nuse; }
	void print () {
	    std::cout << "### " << name << " ###\n";
	    std::cout << "    tot # of builds  : " << nbuild  << "\n"
		      << "    tot # of erasures: " << nerase  << "\n"
		      << "    tot # of uses    : " << nuse    << "\n"
		      << "    max cache size   : " << maxsize << "\n"
		      << "    max # of uses    : " << maxuse
		      << std::endl;
	}
    };
    //
    // Used by a bunch of routines when communicating via MPI.
    //
    struct CopyComTag
    {
        Box box;
        int fabIndex;
        int srcIndex;
	CopyComTag () {}
	explicit CopyComTag (const Box& b, int fidx, int sidx)
	    : box(b), fabIndex(fidx), srcIndex(sidx) {}
	// CopyComTag needs to be sortable if it is used in remote communication.
	// Note that the intersection of two boxes is at most one box.
	// Therefore fabIndex==rhs.fabIndex && srcIndex==rhs.srcIndex means
	// *this == rhs.
	bool operator< (const CopyComTag& rhs) const {
	    return (fabIndex < rhs.fabIndex)
		|| ((fabIndex == rhs.fabIndex) && (srcIndex < rhs.srcIndex));
	}
        //
        // Some typedefs & helper functions used throughout the code.
        //
        typedef std::vector<CopyComTag> CopyComTagsContainer;

        typedef std::map<int,CopyComTagsContainer> MapOfCopyComTagContainers;
    };
    //
    // Some useful typedefs.
    //
    typedef CopyComTag::CopyComTagsContainer CopyComTagsContainer;
    typedef CopyComTag::MapOfCopyComTagContainers MapOfCopyComTagContainers;
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
        // Basic data.
        //
        BoxArray            m_ba;
        DistributionMapping m_dm;
        int                 m_ngrow;
	int                 m_nuse;
        bool                m_cross;
	bool                m_threadsafe_loc;
	bool                m_threadsafe_rcv;
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

    static CacheStats m_FBC_stats;
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
             const DistributionMapping& srcdm,
	     int                        dstng,
	     int                        srcng);

        ~CPC ();

        bool operator== (const CPC& rhs) const;
        bool operator!= (const CPC& rhs) const { return !operator==(rhs); }

        int bytes () const;

        static void FlushCache ();

        BoxArray            m_dstba;
        BoxArray            m_srcba;
        DistributionMapping m_dstdm;
        DistributionMapping m_srcdm;
	int                 m_dstng;
	int                 m_srcng;
        int                 m_nuse;
	bool                m_threadsafe_loc;
	bool                m_threadsafe_rcv;
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

    static CacheStats m_CPC_stats;

    static CPCCacheIter TheCPC (const CPC&          cpc,
                                const FabArrayBase& dst,
                                const FabArrayBase& src);

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
    // Huge box size, larger than any boxes we use
    //
    static IntVect mfiter_huge_box_size;
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
    // Initialize from ParmParse with "fabarray" prefix.
    //
    static void Initialize ();
    static void Finalize ();
    //
    // To maximize thread efficiency we now can decompose things like
    // intersections among boxes into smaller tiles. This sets
    // their maximum size.
    //
    static IntVect comm_tile_size;  // communication tile size

    //
    // The number of FabArrays defined
    //
    static int nFabArrays;
    static int NFabArrays() { return nFabArrays; }
    int FabArrayID() const  { return faID; }

    struct TileArray
    {
	int nuse;
	Array<int> indexMap;	
	Array<int> localIndexMap;
	Array<Box> tileArray;
	TileArray () : nuse(-1) {;}
    };

    const TileArray* getTileArray (const IntVect& tilesize) const;

protected:

    DistributionMapping& ModifyDistributionMap () { return distributionMap; }

    //
    // The data ...
    //
    mutable BoxArray    boxarray;
    DistributionMapping distributionMap;
    Array<int>          indexMap;
    mutable int         n_grow;
    int                 n_comp;
    int                 faID;

    // Key for unique combination of BoxArray and DistributionMapping
    // Note both BoxArray and DistributionMapping are reference counted.
    // Objects with the same references have the same key.
    typedef std::pair<ptrdiff_t,ptrdiff_t> BDKey;

    mutable BDKey m_bdkey;

    BDKey getBDKey () const {
	return std::make_pair(boxarray.getRefID(), distributionMap.getRefID());
    }

    // No operator< for std::pair before C++14
    class CompareBDKey
    {
    public:
	bool operator () (const BDKey& lhs,
			  const BDKey& rhs) const
	    {
		return (lhs.first < rhs.first)
		    || ((lhs.first == rhs.first) && (lhs.second < rhs.second));
	    }
    };

    //
    // Tiling
    //
    // We use tile size as the key for the inner map.
    typedef std::map<IntVect, TileArray, IntVect::Compare> TAMap;
    typedef std::map<BDKey  , TAMap    , CompareBDKey    > TACache;
    //
    static TACache     m_TheTileArrayCache;
    static CacheStats  m_TAC_stats;
    //
    void buildTileArray (const IntVect& tilesize, TileArray& ta) const;
    //
    void flushTileArray (const IntVect& tilesize = IntVect::TheZeroVector()) const;
    // This flushes the entire cache.
    static void flushTileArrayCache ();

    //
    // Keep track of how many FabArrays are built with the same BDKey.
    //
    static std::map<BDKey, int> m_BD_count;
    //
    // clear BD count and TileArray cache associated with this BD, if no other is using this BD. 
    // 
    void clearThisBD ();
    //
    // add the current BD into BD count database
    //
    void addThisBD ();
    //
    struct FabArrayStats
    {
	int  num_fabarrays;
	int  max_num_fabarrays;
	int  max_num_boxarrays;
	int  max_num_ba_use;
	long num_build;
	FabArrayStats () : num_fabarrays(0), max_num_fabarrays(0), max_num_boxarrays(0),
			   max_num_ba_use(1), num_build(0) {;}
	void recordBuild () {
	    ++num_fabarrays;
	    ++num_build;
	    max_num_fabarrays = std::max(max_num_fabarrays, num_fabarrays);
	}
	void recordDelete () {
	    --num_fabarrays;
	}
	void recordMaxNumBoxArrays (int n) {
	    max_num_boxarrays = std::max(max_num_boxarrays, n);
	}
	void recordMaxNumBAUse (int n) {
	    max_num_ba_use = std::max(max_num_ba_use, n);
	}
	void print () {
	    std::cout << "### FabArray ###\n";
	    std::cout << "    tot # of builds       : " << num_build         << "\n"
		      << "    max # of FabArrays    : " << max_num_fabarrays << "\n"
		      << "    max # of BoxArrays    : " << max_num_boxarrays << "\n"
		      << "    max # of BoxArray uses: " << max_num_ba_use
		      << std::endl;
	}
    };
    static FabArrayStats m_FA_stats;
};

class MFIter
{
public:
    enum Flags {
        Tiling        = 0x01,
        NoSharing     = 0x02, // For OMP only. If on, all threads loops over all grids.
        NoTeamBarrier = 0x04  // For Team only. If on, there is no barrier in MFIter dtor.
    };  // All these flags are off by default.
    //
    // Construct a MFIter.
    //
    // The default is no tiling
    explicit MFIter (const FabArrayBase& fabarray,
		     unsigned char       flags_=0);
    // tiling w/ default size, IntVect FabArrayBase::mfiter_tile_size
    explicit MFIter (const FabArrayBase& fabarray, 
		     bool                do_tiling); 
    // tiling with explicit size and flags
    explicit MFIter (const FabArrayBase& fabarray, 
		     const IntVect&      tilesize,
		     unsigned char       flags_=0);
    // tiling with explicit tilearray and flags
    explicit MFIter (const FabArrayBase&            fabarray, 
		     const FabArrayBase::TileArray* pta_,
		     unsigned char                  flags_=0);
    // dtor
    ~MFIter ();
    //
    // Returns the tile Box at the current index.
    //
    Box tilebox () const;
    //
    // Returns the dir-nodal (or all nodal if dir<0) Box at the current index.
    //
    Box nodaltilebox (int dir=-1) const;
    //
    // Returns the tile box at the current index grown to include ghost cells.
    //
    Box growntilebox (int ng=-1000000) const;
    //
    //  Returns the dir-nodal (or all nodal if dir<0) box grown to include ghost cells.
    //
    Box grownnodaltilebox (int dir=-1, int ng=-1000000) const;
    //
    // Returns the valid Box that current tile resides.
    //
    Box validbox () const { return fabArray.box(pta->indexMap[currentIndex]); }
    //
    // Returns the Box of the FAB at which we currently point.
    //
    Box fabbox () const { return fabArray.fabbox(pta->indexMap[currentIndex]); }
    //
    // Increments iterator to the next tile we own.
    //
    void operator++ () { ++currentIndex;}
    //
    // Is the iterator valid i.e. is it associated with a FAB?
    //
    bool isValid () { return currentIndex < endIndex; }
    //
    // The index into the underlying BoxArray of the current FAB.
    //
    int index () const { return pta->indexMap[currentIndex]; }
    //
    // local index into the vector of fab pointers, m_fabs_v
    //
    int LocalIndex () const { return pta->localIndexMap[currentIndex]; }
    //
    // Constant reference to FabArray over which we're iterating.
    //
    const FabArrayBase& theFabArrayBase () const { return fabArray; }

private:

    const FabArrayBase& fabArray;
    const FabArrayBase::TileArray* pta;

    unsigned char flags;
    int           currentIndex;
    int           beginIndex;
    int           endIndex;
    IndexType     typ;

    void Initialize ();
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
              FabAlloc        mem_mode = Fab_allocate,
	      const IntVect&  nodal = IntVect::TheZeroVector());

    FabArray (const BoxArray&            bxs,
              int                        nvar,
              int                        ngrow,
              const DistributionMapping& dm,
              FabAlloc                   mem_mode = Fab_allocate,
	      const IntVect&             nodal = IntVect::TheZeroVector());
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
                 FabAlloc        mem_mode,
		 const IntVect&  nodal = IntVect::TheZeroVector());

    void define (const BoxArray&            bxs,
                 int                        nvar,
                 int                        ngrow,
		 const DistributionMapping& dm,
                 FabAlloc                   mem_mode,
		 const IntVect&             nodal = IntVect::TheZeroVector());
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
    // is restricted to the valid region of destination and the 
    // valid+src_nghost region of source.
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
    void copy (const FabArray<FAB>& src,
               int                  src_comp,
               int                  dest_comp,
               int                  num_comp,
	       int                  src_nghost,
	       int                  dst_nghost,
               CpOp                 op = FabArrayBase::COPY);
    //
    // Copies the values contained in the intersection of the
    // valid + nghost region of this FabArray with the FAB dest into dest.
    //
    void copy (FAB& dest,
	       int  nghost = 0) const;
    //
    // Copies the values contained in the intersection of the
    // valid + nghost region of this FabArray with the FAB dest and the Box
    // subbox into that subregion of dest.
    //
    void copy (FAB&       dest,
               const Box& subbox,
	       int        nghost = 0) const;
    //
    // Copies the values contained in the intersection of the
    // num_comp component valid + nghost region of this FabArray, starting at
    // component src_comp, with the FAB dest into dest, starting at
    // component dest_comp in dest.
    //
    void copy (FAB& dest,
               int  src_comp,
               int  dest_comp,
               int  num_comp,
	       int  nghost = 0) const;
    //
    // Copies the values contained in the intersection of the
    // num_comp component valid + nghost region of this FabArray, starting at
    // component src_comp, with the FAB dest and the Box subbox, into
    // dest, starting at component dest_comp in dest.
    //
    void copy (FAB&       dest,
               const Box& subbox,
               int        src_comp,
               int        dest_comp,
               int        num_comp,
	       int        nghost = 0) const;

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

    //
    // Move FABs in this FabArray to different MPI ranks.
    //

    struct FABMoves {
      int distMapIndex, fromRank, toRank, seqNum;
    };

    int MoveFabs (const Array<int> &newDistMapArray);
    static void MoveAllFabs (const Array<int> &newDistMapArray);

protected:
    //
    // Helper function for define().
    //
    void defineDoit (const BoxArray&            bxs,
                     int                        nvar,
                     int                        ngrow,
                     FabAlloc                   mem_mode,
                     const DistributionMapping* dm,
		     const IntVect&             nodal);
    //
    // The data.
    //
    std::vector<FAB*> m_fabs_v;

    static std::map<int, std::map<int, FabArray<FAB> *> > allocatedFAPointers;
           // <ngrids to find distmap, <FabArrayID, FabArray *> >

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

template<typename T>
void
FabArrayBase::PostRcvs (const std::map<int,int>&               m_RcvVols,
                        T*&                                    the_recv_data,
                        Array<T*>&                             recv_data,
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

    BL_ASSERT((TotalRcvsVolume*sizeof(T)) < std::numeric_limits<int>::max());

    the_recv_data = static_cast<T*>(BoxLib::The_Arena()->alloc(TotalRcvsVolume*sizeof(T)));

    int Offset = 0;

    for (std::map<int,int>::const_iterator it = m_RcvVols.begin(),
             End = m_RcvVols.end();
         it != End;
         ++it)
    {
        const int N = it->second*ncomp;

        BL_ASSERT(N < std::numeric_limits<int>::max());

        recv_data.push_back(&the_recv_data[Offset]);
        recv_from.push_back(it->first);
        recv_reqs.push_back(ParallelDescriptor::Arecv(recv_data.back(),N,it->first,SeqNum).req());

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
    BL_COMM_PROFILE_WAITSOME(BLProfiler::Waitall, send_reqs, N_snds, indx, stats, false);

    BL_MPI_REQUIRE( MPI_Waitall(N_snds, send_reqs.dataPtr(), stats.dataPtr()) );

    BL_COMM_PROFILE_WAITSOME(BLProfiler::Waitall, send_reqs, N_snds, indx, stats, false);

    for (int i = 0; i < N_snds; i++)
        BoxLib::The_Arena()->free(send_data[i]);
#endif /*BL_USE_MPI*/
}

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
    if (li < m_fabs_v.size() && m_fabs_v[li] != 0) {
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
    BL_ASSERT(mfi.LocalIndex() < indexMap.size());
    return *m_fabs_v[mfi.LocalIndex()];
}

template <class FAB>
FAB&
FabArray<FAB>::operator[] (const MFIter& mfi)
{
    BL_ASSERT(mfi.LocalIndex() < indexMap.size());
    return *m_fabs_v[mfi.LocalIndex()];
}

template <class FAB>
const FAB&
FabArray<FAB>::operator[] (int K) const
{
    int li = localindex(K);
    BL_ASSERT(li >=0 && li < indexMap.size());
    return *m_fabs_v[li];
}

template <class FAB>
FAB&
FabArray<FAB>::operator[] (int K)
{
    int li = localindex(K);
    BL_ASSERT(li >=0 && li < indexMap.size());
    return *m_fabs_v[li];
}

template <class FAB>
void
FabArray<FAB>::clear ()
{
    clearThisBD();

    typename std::map<int, std::map<int, FabArray<FAB> *> >::iterator afapIter = 
                 FabArray<FAB>::allocatedFAPointers.find(distributionMap.size());
    if(afapIter == FabArray<FAB>::allocatedFAPointers.end()) {
#ifdef DEBUG_AFAP
      if(ParallelDescriptor::IOProcessor() && distributionMap.size() > 0) {
        std::cout << "**** In FabArray::clear():: map not found:  size = "
                  << distributionMap.size() << std::endl;
      }
#endif
    } else {
      std::map<int, FabArray<FAB> *> &faPtrCachedMap = afapIter->second;
      faPtrCachedMap.erase(faID);
    }

    for(Iterator it = m_fabs_v.begin(); it != m_fabs_v.end(); ++it) {
	delete *it;
    }
    
    m_fabs_v.clear();
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
{
    m_FA_stats.recordBuild();
}

template <class FAB>
FabArray<FAB>::FabArray (const BoxArray& bxs,
                         int             nvar,
                         int             ngrow,
                         FabAlloc        alloc,
			 const IntVect&  nodal)
{
    m_FA_stats.recordBuild();
    define(bxs,nvar,ngrow,alloc,nodal);
}

template <class FAB>
FabArray<FAB>::FabArray (const BoxArray&            bxs,
                         int                        nvar,
                         int                        ngrow,
                         const DistributionMapping& dm,
                         FabAlloc                   alloc,
			 const IntVect&             nodal)
{
    m_FA_stats.recordBuild();
    define(bxs,nvar,ngrow,dm,alloc,nodal);
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
                           const DistributionMapping* dm,
			   const IntVect&             nodal)
{
    BL_ASSERT(ngrow >= 0);
    BL_ASSERT(boxarray.size() == 0);

    n_grow = ngrow;
    n_comp = nvar;

    boxarray.define(bxs);

    if (nodal != IntVect::TheZeroVector()) {
	if (nodal == IntVect::TheUnitVector()) {
	    boxarray.surroundingNodes();
	} else {
	    for (int i=0; i<BL_SPACEDIM; ++i) {
		if (nodal[i]) {
		    boxarray.surroundingNodes(i);
		} 
	    }
	}
    }

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

    for(int i = 0, N = boxarray.size(); i < N; ++i) {
      if(distributionMap[i] == MyProc) {
        indexMap.push_back(i);
      }
    }

    addThisBD();
    
    if(alloc == Fab_allocate) {
      AllocFabs();

      typename std::map<int, std::map<int, FabArray<FAB> *> >::iterator afapIter = 
                   FabArray<FAB>::allocatedFAPointers.find(distributionMap.size());
      if(afapIter == FabArray<FAB>::allocatedFAPointers.end()) {
	int dmapSize(distributionMap.size());
	std::map<int, FabArray<FAB> *> tempMap;
	tempMap.insert(std::make_pair(faID, this));
        FabArray<FAB>::allocatedFAPointers.insert(std::make_pair(dmapSize, tempMap));

      } else {
        std::map<int, FabArray<FAB> *> &faPtrCachedMap = afapIter->second;
        faPtrCachedMap.insert(std::make_pair(faID, this));
      }
    }
}

template <class FAB>
void
FabArray<FAB>::define (const BoxArray& bxs,
                       int             nvar,
                       int             ngrow,
                       FabAlloc        alloc,
		       const IntVect&  nodal)
{
    defineDoit(bxs,nvar,ngrow,alloc,0,nodal);
}

template <class FAB>
void
FabArray<FAB>::define (const BoxArray&            bxs,
                       int                        nvar,
                       int                        ngrow,
                       const DistributionMapping& dm,
                       FabAlloc                   alloc,
		       const IntVect&             nodal)
{
    defineDoit(bxs,nvar,ngrow,alloc,&dm,nodal);
}

template <class FAB>
void
FabArray<FAB>::AllocFabs ()
{
    m_fabs_v.reserve(indexMap.size());

    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        const Box& tmp = BoxLib::grow(fai.validbox(), n_grow);
	m_fabs_v.push_back(new FAB(tmp, n_comp));
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

    if (m_fabs_v.size() == 0) {
	m_fabs_v.resize(indexMap.size());
    }

    m_fabs_v[localindex(boxno)] = elem;
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

    if (m_fabs_v.size() == 0) {
	m_fabs_v.resize(indexMap.size());
    }

    m_fabs_v[mfi.LocalIndex()] = elem;
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
void
FabArray<FAB>::copy (const FabArray<FAB>& src,
                     int                  scomp,
                     int                  dcomp,
                     int                  ncomp,
		     int                  snghost,
		     int                  dnghost,
                     CpOp                 op)
{
    BL_PROFILE("FabArray::copy()");

    if (size() == 0 || src.size() == 0) return;

    BL_ASSERT(op == FabArrayBase::COPY || op == FabArrayBase::ADD);
    BL_ASSERT(boxArray().ixType() == src.boxArray().ixType());

    BL_ASSERT(src.nGrow() >= snghost);
    BL_ASSERT(    nGrow() >= dnghost);

    if ((src.boxArray().ixType().cellCentered() || op == FabArrayBase::COPY) &&
        (boxarray == src.boxarray && distributionMap == src.distributionMap))
    {
        //
        // Short-circuit full intersection code if we're doing copy()s or if
        // we're doing plus()s on cell-centered data.  Don't do plus()s on
        // non-cell-centered data this simplistic way.
        //
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter fai(*this,true); fai.isValid(); ++fai)
        {
            const Box& bx = fai.tilebox();

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

    const CPC cpc(boxarray, src.boxarray, distributionMap, src.distributionMap, dnghost, snghost);

    FabArrayBase::CPCCacheIter cache_it = FabArrayBase::TheCPC(cpc, *this, src);

    BL_ASSERT(cache_it != FabArrayBase::m_TheCopyCache.end());

    const CPC& thecpc = cache_it->second;

    if (ParallelDescriptor::NProcs() == 1)
    {
        //
        // There can only be local work to do.
        //
	int N_loc = (*thecpc.m_LocTags).size();
#ifdef _OPENMP
#pragma omp parallel for if (thecpc.m_threadsafe_loc)
#endif
	for (int i=0; i<N_loc; ++i)
        {
            const CopyComTag& tag = (*thecpc.m_LocTags)[i];

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

        Array<MPI_Status>  stats;
        Array<int>         recv_from;
        Array<value_type*> recv_data;
        Array<MPI_Request> recv_reqs;
        //
        // Post rcvs. Allocate one chunk of space to hold'm all.
        //
        value_type* the_recv_data = 0;

        FabArrayBase::PostRcvs(*thecpc.m_RcvVols,the_recv_data,recv_data,recv_from,recv_reqs,NC,SeqNum);

	//
	// Post send's
	// 
	const int N_snds = thecpc.m_SndTags->size();

	Array<value_type*>                 send_data;
	Array<int>                         send_N;
	Array<int>                         send_rank;
	Array<const CopyComTagsContainer*> send_cctc;

	send_data.reserve(N_snds);
	send_N   .reserve(N_snds);
	send_rank.reserve(N_snds);
	send_cctc.reserve(N_snds);

        for (MapOfCopyComTagContainers::const_iterator m_it = thecpc.m_SndTags->begin(),
                 m_End = thecpc.m_SndTags->end();
             m_it != m_End;
             ++m_it)
        {
            std::map<int,int>::const_iterator vol_it = thecpc.m_SndVols->find(m_it->first);

            BL_ASSERT(vol_it != thecpc.m_SndVols->end());

            const int N = vol_it->second*NC;

            BL_ASSERT(N < std::numeric_limits<int>::max());

            value_type* data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(N*sizeof(value_type)));
 
	    send_data.push_back(data);
	    send_N   .push_back(N);
	    send_rank.push_back(m_it->first);
	    send_cctc.push_back(&(m_it->second));
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int j=0; j<N_snds; ++j)
	{
	    value_type* dptr = send_data[j];
	    BL_ASSERT(dptr != 0);

	    const CopyComTagsContainer& cctc = *send_cctc[j];

            for (CopyComTagsContainer::const_iterator it = cctc.begin();
                 it != cctc.end(); ++it)
            {
                const Box& bx = it->box;
                src[it->srcIndex].copyToMem(bx,SC,NC,dptr);
                const int Cnt = bx.numPts()*NC;
                dptr += Cnt;
            }
	}

	Array<MPI_Request> send_reqs;

	if (FabArrayBase::do_async_sends)
	{
	    send_reqs.reserve(N_snds);
	    for (int j=0; j<N_snds; ++j)
	    {
                send_reqs.push_back(ParallelDescriptor::Asend
				    (send_data[j],send_N[j],send_rank[j],SeqNum).req());
            }
	} else {
	    for (int j=0; j<N_snds; ++j)
	    {
                ParallelDescriptor::Send(send_data[j],send_N[j],send_rank[j],SeqNum);
                BoxLib::The_Arena()->free(send_data[j]);
            }
        }

        //
        // Do the local work.  Hope for a bit of communication/computation overlap.
        //
	int N_loc = (*thecpc.m_LocTags).size();
#ifdef _OPENMP
#pragma omp parallel for if (thecpc.m_threadsafe_loc)
#endif
	for (int j=0; j<N_loc; ++j)
        {
            const CopyComTag& tag = (*thecpc.m_LocTags)[j];

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
	    Array<const CopyComTagsContainer*> recv_cctc;
	    recv_cctc.reserve(N_rcvs);

	    for (int k = 0; k < N_rcvs; k++)
	    {
		MapOfCopyComTagContainers::const_iterator m_it = thecpc.m_RcvTags->find(recv_from[k]);
		BL_ASSERT(m_it != thecpc.m_RcvTags->end());
	    
		recv_cctc.push_back(&(m_it->second));
	    }

	    stats.resize(N_rcvs);
	    BL_MPI_REQUIRE( MPI_Waitall(N_rcvs, recv_reqs.dataPtr(), stats.dataPtr()) );
	    
#ifdef _OPENMP
#pragma omp parallel if (thecpc.m_threadsafe_rcv)
#endif
        {
	    FAB fab;

#ifdef _OPENMP
#pragma omp for
#endif
	    for (int k = 0; k < N_rcvs; k++)
	    {
		const value_type* dptr = recv_data[k];
		BL_ASSERT(dptr != 0);
		
		const CopyComTagsContainer& cctc = *recv_cctc[k];
		
		for (CopyComTagsContainer::const_iterator it = cctc.begin();
		     it != cctc.end(); ++it)
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
	}
	
        BoxLib::The_Arena()->free(the_recv_data);
	
        if (FabArrayBase::do_async_sends && !thecpc.m_SndTags->empty())
            FabArrayBase::GrokAsyncSends(thecpc.m_SndTags->size(),send_reqs,send_data,stats);

        ipass     += NC;
        SC        += NC;
        DC        += NC;
        NCompLeft -= NC;
    }

    return;

#endif /*BL_USE_MPI*/
}

template <class FAB>
void
FabArray<FAB>::copy (const FabArray<FAB>& src,
                     int                  scomp,
                     int                  dcomp,
                     int                  ncomp,
                     CpOp                 op)
{
    copy(src,scomp,dcomp,ncomp,0,0,op);
}

template <class FAB>
void
FabArray<FAB>::copy (const FabArray<FAB>& src, CpOp op)
{
    copy(src,0,0,nComp(),0,0,op);
}

//
// Copies to FABs, note that destination is first arg.
//

template <class FAB>
void
FabArray<FAB>::copy (FAB& dest,
		     int  nghsot) const
{
    copy(dest, dest.box(), 0, 0, dest.nComp(), nghsot);
}

template <class FAB>
void
FabArray<FAB>::copy (FAB&       dest,
                     const Box& subbox,
		     int        nghost) const
{
    copy(dest, subbox, 0, 0, dest.nComp(), nghost);
}

template <class FAB>
void
FabArray<FAB>::copy (FAB& dest,
                     int  scomp,
                     int  dcomp,
                     int  ncomp,
		     int  nghost) const
{
    copy(dest, dest.box(), scomp, dcomp, ncomp, nghost);
}

template <class FAB>
void
FabArray<FAB>::copy (FAB&       dest,
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
	    const Box& bx = BoxLib::grow(boxarray[j],nghost);
	    const Box& destbox = bx & subbox;
	    if (destbox.ok())
            {
                dest.copy(get(j),destbox,scomp,destbox,dcomp,ncomp);
            }
        }

        return;
    }

    FAB ovlp;

    for (int i = 0, N = size(); i < N; i++)
    {
	const Box& srcbox = BoxLib::grow(boxarray[i],nghost);
	const Box& bx = srcbox & subbox;
	if (bx.ok())
        {
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
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter fai(*this,true); fai.isValid(); ++fai)
    {
	Box bx = fai.growntilebox();
        get(fai).setVal(val, bx, 0, n_comp);
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

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter fai(*this,true); fai.isValid(); ++fai)
    {
	Box bx = fai.growntilebox(nghost);
        get(fai).setVal(val, bx, comp, ncomp);
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
FabArray<FAB>::MoveAllFabs (const Array<int> &newDistMapArray)
{
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "FabArray<FAB>::MoveAllFabs:  " << allocatedFAPointers.size()
              << " cached DistributionMap(s)." << std::endl;
  }
  FabArray<FAB> *lastFAPtr = 0;

  typename std::map<int, std::map<int, FabArray<FAB> *> >::iterator afapIter = 
               FabArray<FAB>::allocatedFAPointers.find(newDistMapArray.size());
  if(afapIter == FabArray<FAB>::allocatedFAPointers.end()) {
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "FabArray<FAB>::MoveAllFabs:  no allocated pointers for new "
                << "distribution map array." << std::endl;
    }
  } else {
    int nFabsMoved(0);
    std::map<int, FabArray<FAB> *> &faPtrCachedMap = afapIter->second;
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "FabArray<FAB>::MoveAllFabs:  moving FABs for "
                << faPtrCachedMap.size() << " FabArray(s)." << std::endl;
    }

    for(typename std::map<int, FabArray<FAB> *>::iterator it = faPtrCachedMap.begin();
        it != faPtrCachedMap.end(); ++it)
    {
      if(it->second->ok() == false) {
        BoxLib::Abort("it not ok");
      }
      nFabsMoved = it->second->MoveFabs(newDistMapArray);  // just keep the last one
      if( ! it->second->ok()) {
        BoxLib::Abort("_here 00:  it not ok");
      }
      lastFAPtr = it->second;
    }
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "FabArray<FAB>::MoveAllFabs:  moved " << nFabsMoved
	        << " FAB(s) for each FabArray." << std::endl;
    }

    // ---- flush caches
    DistributionMapping::FlushCache();
    FabArrayBase::FlushSICache();
    FabArrayBase::CPC::FlushCache();
    //Geometry::FlushPIRMCache();  // ---- needs to be done externally

    if(lastFAPtr != 0) {
      lastFAPtr->ModifyDistributionMap().ReplaceCachedProcessorMap(newDistMapArray);
    }
  }
}


template <class FAB>
int
FabArray<FAB>::MoveFabs (const Array<int> &newDistMapArray)
{
#if BL_USE_MPI
  BL_PROFILE("FabArray<FAB>::MoveFabs()");

  int myProc(ParallelDescriptor::MyProc());

  // ---- check validity of newDistMapArray
  if(newDistMapArray.size() != distributionMap.size()) {
    std::cout << "ndma.size  dm.size = " << newDistMapArray.size() << "  "
              << distributionMap.size() << std::endl;
    BoxLib::Abort("**** Error:  bad newDistMap:0");
  }
  if(newDistMapArray[newDistMapArray.size() - 1] != myProc) {
    BoxLib::Abort("**** Error:  bad newDistMap:1");
  }
  for(int idm(0); idm < newDistMapArray.size(); ++idm) {
    if(newDistMapArray[idm] < 0 ||
       newDistMapArray[idm] > ParallelDescriptor::NProcs() - 1)
    {
      BoxLib::Abort("**** Error:  bad newDistMap:2");
    }
  }
  if(newDistMapArray == distributionMap.ProcessorMap()) {
    return 0;
  }

  // ---- determine which fabs to move
  std::vector<FABMoves> fabMoves;
  for(int iM(0); iM < distributionMap.size() - 1; ++iM) {  // ---- -1 skips the sentinel
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
    if(indexMap.size() != m_fabs_v.size()) {
      BoxLib::Abort("**** Error:  indexMap size != m_fabs_v.size()");
    }
    for(int iim(0); iim < indexMap.size(); ++iim) {
      tIFiter = tempIndexFABs.find(indexMap[iim]);
      if(tIFiter == tempIndexFABs.end()) {
	tempIndexFABs.insert(std::pair<int, FAB *>(indexMap[iim], m_fabs_v[iim]));
      } else {
        BoxLib::Abort("**** Error:  index not in indexMap.");
      }
    }
  }

  // ---- move the fabs
  Array<MPI_Request> recvReqs, sendReqs;
  int nFabsSent(0);

  for(int imoves(0); imoves < fabMoves.size(); ++imoves) {
    FABMoves &moveThisFab = fabMoves[imoves];
    int dmi(moveThisFab.distMapIndex);

    if(myProc == moveThisFab.toRank) {   // ---- receive fab(s)
      const Box &tmpbox = BoxLib::grow(boxarray[dmi], n_grow);
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
        BoxLib::Abort("**** Error:  index not in tempIndexFABs.");
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

  // ---- we could defer this, not do it, or just do it the first time through
  ParallelDescriptor::ReduceIntSum(nFabsSent);


  // ---- wait for all the data to move
  // ---- we could defer Waitall for multiple calls to MoveFabs with a multi-step process
  Array<MPI_Status>  recvStats(recvReqs.size()), sendStats(sendReqs.size());

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
        BoxLib::Abort("**** Error:  index not in tempIndexFABs when deleting.");
      } else {
	delete tIFiter->second;
	tempIndexFABs.erase(tIFiter);
      }
    }
  }

  // ---- reconstruct the index and fab vectors
  indexMap.clear();
  m_fabs_v.clear();
  for(tIFiter = tempIndexFABs.begin(); tIFiter != tempIndexFABs.end(); ++ tIFiter) {
    indexMap.push_back(tIFiter->first);
    m_fabs_v.push_back(tIFiter->second);
  }

  return nFabsSent;
#else
  return 0;
#endif
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
	    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int), ParallelDescriptor::MyProc(),
	                    BLProfiler::BeforeCall());

            BL_MPI_REQUIRE( MPI_Alltoall(RcvsArray.dataPtr(),
                                         1,
                                         ParallelDescriptor::Mpi_typemap<int>::type(),
                                         SndsArray.dataPtr(),
                                         1,
                                         ParallelDescriptor::Mpi_typemap<int>::type(),
                                         ParallelDescriptor::Communicator()) );

	    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int), ParallelDescriptor::MyProc(),
	                    BLProfiler::AfterCall());

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
	    BL_COMM_PROFILE(BLProfiler::Alltoallv, nrcvs * sizeof(int), ParallelDescriptor::MyProc(),
	                    BLProfiler::BeforeCall());

            BL_MPI_REQUIRE( MPI_Alltoallv(cd_that_i_need.dataPtr(),
                                          scnts.dataPtr(),
                                          sdispls.dataPtr(),
                                          ParallelDescriptor::Mpi_typemap<int>::type(),
                                          cd_others_need.dataPtr(),
                                          rcnts.dataPtr(),
                                          rdispls.dataPtr(),
                                          ParallelDescriptor::Mpi_typemap<int>::type(),
                                          ParallelDescriptor::Communicator()) );

	    BL_COMM_PROFILE(BLProfiler::Alltoallv, nrcvs * sizeof(int), ParallelDescriptor::MyProc(),
	                    BLProfiler::AfterCall());
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

        BL_COMM_PROFILE(BLProfiler::Waitall, sizeof(value_type), BLProfiler::BeforeCall(), N_snds);
        BL_MPI_REQUIRE( MPI_Waitall(N_snds, send_reqs.dataPtr(), stats.dataPtr()) );
        BL_COMM_PROFILE(BLProfiler::Waitall, sizeof(value_type), BLProfiler::AfterCall(), N_snds);

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
	int N_loc = (*TheSI.m_LocTags).size();
#ifdef _OPENMP
#pragma omp parallel for if (TheSI.m_threadsafe_loc)
#endif
	for (int i=0; i<N_loc; ++i)
        {
            const CopyComTag& tag = (*TheSI.m_LocTags)[i];

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
    Array<value_type*> recv_data;
    Array<MPI_Request> recv_reqs;
    //
    // Post rcvs. Allocate one chunk of space to hold'm all.
    //
    value_type* the_recv_data = 0;

    FabArrayBase::PostRcvs(*TheSI.m_RcvVols,the_recv_data,recv_data,recv_from,recv_reqs,ncomp,SeqNum);

    //
    // Post send's
    //
    const int N_snds = TheSI.m_SndTags->size();

    Array<value_type*>                 send_data;
    Array<int>                         send_N;
    Array<int>                         send_rank;
    Array<const CopyComTagsContainer*> send_cctc;
    
    send_data.reserve(N_snds);
    send_N   .reserve(N_snds);
    send_rank.reserve(N_snds);
    send_cctc.reserve(N_snds);

    for (MapOfCopyComTagContainers::const_iterator m_it = TheSI.m_SndTags->begin(),
             m_End = TheSI.m_SndTags->end();
         m_it != m_End;
         ++m_it)
    {
	std::map<int,int>::const_iterator vol_it = TheSI.m_SndVols->find(m_it->first);

        BL_ASSERT(vol_it != TheSI.m_SndVols->end());

        const int N = vol_it->second*ncomp;

        BL_ASSERT(N < std::numeric_limits<int>::max());

        value_type* data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(N*sizeof(value_type)));

	send_data.push_back(data);
	send_N   .push_back(N);
	send_rank.push_back(m_it->first);
	send_cctc.push_back(&(m_it->second));
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=0; i<N_snds; ++i)
    {
	value_type* dptr = send_data[i];
	BL_ASSERT(dptr != 0);

	const CopyComTagsContainer& cctc = *send_cctc[i];

	for (CopyComTagsContainer::const_iterator it = cctc.begin();
		 it != cctc.end(); ++it)
        {
            BL_ASSERT(distributionMap[it->srcIndex] == ParallelDescriptor::MyProc());
            const Box& bx = it->box;
            get(it->srcIndex).copyToMem(bx,scomp,ncomp,dptr);
            const int Cnt = bx.numPts()*ncomp;
            dptr += Cnt;
        }
    }

    Array<MPI_Request> send_reqs;

    if (FabArrayBase::do_async_sends)
    {
	send_reqs.reserve(N_snds);
	for (int i=0; i<N_snds; ++i) {
	    send_reqs.push_back(ParallelDescriptor::Asend
				(send_data[i],send_N[i],send_rank[i],SeqNum).req());
	}
    } else {
	for (int i=0; i<N_snds; ++i) {
	    ParallelDescriptor::Send(send_data[i],send_N[i],send_rank[i],SeqNum);
	    BoxLib::The_Arena()->free(send_data[i]);
	}
    }

    //
    // Do the local work.  Hope for a bit of communication/computation overlap.
    //
    int N_loc = (*TheSI.m_LocTags).size();
#ifdef _OPENMP
#pragma omp parallel for if (TheSI.m_threadsafe_loc)
#endif
    for (int i=0; i<N_loc; ++i)
    {
        const CopyComTag& tag = (*TheSI.m_LocTags)[i];

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
	Array<const CopyComTagsContainer*> recv_cctc;
	recv_cctc.reserve(N_rcvs);

	for (int k = 0; k < N_rcvs; k++) 
	{
	    MapOfCopyComTagContainers::const_iterator m_it = TheSI.m_RcvTags->find(recv_from[k]);
	    BL_ASSERT(m_it != TheSI.m_RcvTags->end());
	    
	    recv_cctc.push_back(&(m_it->second));
	}	

	stats.resize(N_rcvs);
	BL_MPI_REQUIRE( MPI_Waitall(N_rcvs, recv_reqs.dataPtr(), stats.dataPtr()) );

#ifdef _OPENMP
#pragma omp parallel for if (TheSI.m_threadsafe_rcv)
#endif
	for (int k = 0; k < N_rcvs; k++) 
	{
	    value_type*  dptr = recv_data[k];
	    BL_ASSERT(dptr != 0);

	    const CopyComTagsContainer& cctc = *recv_cctc[k];

	    for (CopyComTagsContainer::const_iterator it = cctc.begin();
		 it != cctc.end(); ++it)
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


template <typename FAB> std::map<int, std::map<int, FabArray<FAB> *> > FabArray<FAB>::allocatedFAPointers;

#endif /*BL_FABARRAY_H*/
