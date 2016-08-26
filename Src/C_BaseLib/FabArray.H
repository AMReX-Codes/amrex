
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
#include <BLPgas.H>
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
#include <Utility.H>
#include <ccse-mpi.H>
#include <BLProfiler.H>
#include <Periodicity.H>

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

class MFIter;
class MFGhostIter;

class FabArrayBase
{
    friend class MFIter;
    friend class MFGhostIter;

public:

    FabArrayBase ();

    virtual ~FabArrayBase();
    //
    // Returns the grow factor that defines the region of definition.
    //
    int nGrow () const { return n_grow; }
    //
    // Returns number of variables associated with each point (nvar).
    //
    int nComp () const { return n_comp; }
    //
    ParallelDescriptor::Color color () const { return distributionMap.color(); }
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
    int local_size () const { return indexArray.size(); }
    //
    // Returns constant reference to indices in the FabArray that we have access
    //
    const Array<int> &IndexArray () const { return indexArray; }
    //
    // Returns local index in the vector of FABs.
    //
    int localindex (int K) const { 
	std::vector<int>::const_iterator low
	    = std::lower_bound(indexArray.begin(), indexArray.end(), K);
	if (low != indexArray.end() && *low == K) {
	    return low - indexArray.begin();
	}
	else {
	    return -1;
	}
    }
    //
    // Returns constant reference to associated DistributionMapping.
    //
    const DistributionMapping& DistributionMap () const { return distributionMap; }
    //
    // Some static member templates used throughout the code.
    //
    template<typename T>
    static void WaitForAsyncSends (int                 N_snds,
                                   Array<MPI_Request>& send_reqs,
                                   Array<T*>&          send_data,
                                   Array<MPI_Status>&  stats);

#ifdef BL_USE_UPCXX
    template<typename T>
    static void WaitForAsyncSends_PGAS (int                 N_snds,
                                        Array<T*>&          send_data,
                                        upcxx::event*       send_event,
                                        volatile int*       send_counter);
#endif

    template<typename T>
    static void PostRcvs (const std::map<int,int>&               m_RcvVols,
                          T*&                                    the_recv_data,
                          Array<T*>&                             recv_data,
                          Array<int>&                            recv_from,
                          Array<MPI_Request>&                    recv_reqs,
                          int                                    ncomp,
                          int                                    SeqNum,
			  MPI_Comm   comm = ParallelDescriptor::Communicator());

#ifdef BL_USE_MPI3
    template<typename T>
    static void PostRcvs_MPI_Onesided (
			  const std::map<int,int>&               m_RcvVols,
                          T*&                                    the_recv_data,
                          Array<T*>&                             recv_data,
                          Array<int>&                            recv_from,
                          Array<MPI_Request>&                    recv_reqs,
                          Array<MPI_Aint>          &             recv_disp,
                          int                                    ncomp,
                          int                                    SeqNum,
			  MPI_Win &                              win);
#endif
  
#ifdef BL_USE_UPCXX
    template<typename T>
    static void PostRcvs_PGAS (const std::map<int,int>&               m_RcvVols,
                               T*&                                    the_recv_data,
                               Array<T*>&                             recv_data,
                               Array<int>&                            recv_from,
                               int                                    ncomp,
                               int                                    SeqNum,
                               upcxx::event*                          recv_event);
#endif

    //
    struct CacheStats
    {
	int         size;     // current size: nbuild - nerase
	int         maxsize;  // highest water mark of size
	int         maxuse;   // max # of uses of a cached item
	long        nuse;     // # of uses of the whole cache
	long        nbuild;   // # of build operations
	long        nerase;   // # of erase operations
	long        bytes;
	long        bytes_hwm;
	std::string name;     // name of the cache
	CacheStats (const std::string& name_) 
	    : size(0),maxsize(0),maxuse(0),nuse(0),nbuild(0),nerase(0),
	      bytes(0L),bytes_hwm(0L),name(name_) {;}
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
	Box dbox;
        Box sbox;
        int dstIndex;
        int srcIndex;
	CopyComTag () {}
	CopyComTag (const Box& db, const Box& sb, int didx, int sidx)
	    : dbox(db), sbox(sb), dstIndex(didx), srcIndex(sidx) {}
	bool operator< (const CopyComTag& rhs) const {
	    return (srcIndex < rhs.srcIndex) || ((srcIndex == rhs.srcIndex) && (
                   (dstIndex < rhs.dstIndex) || ((dstIndex == rhs.dstIndex) && (
                   (IntVect::Compare()(dbox.smallEnd(),rhs.dbox.smallEnd())
		               || ((dbox.smallEnd() == rhs.dbox.smallEnd()) && (
                   (IntVect::Compare()(sbox.smallEnd(),rhs.sbox.smallEnd())))))))));
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
    static long bytesOfMapOfCopyComTagContainers (const MapOfCopyComTagContainers&);

    // Key for unique combination of BoxArray and DistributionMapping
    // Note both BoxArray and DistributionMapping are reference counted.
    // Objects with the same references have the same key.
    typedef std::pair<ptrdiff_t,ptrdiff_t> BDKey;

    BDKey getBDKey () const {
	return std::make_pair(boxarray.getRefID(), distributionMap.getRefID());
    }

    void updateBDKey ();

    //
    // Tiling
    //
    struct TileArray
    {
	int nuse;
	Array<int> indexMap;	
	Array<int> localIndexMap;
	Array<Box> tileArray;
	TileArray () : nuse(-1) {;}
	long bytes () const;
    };

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
    // Default tilesize in MFIter
    //
    static IntVect mfiter_tile_size;
    //
    // Default tilesize in MFGhostIter
    //
    static IntVect mfghostiter_tile_size;
    //
    // The maximum number of components to copy() at a time.
    //
    static int MaxComp;
    //
    // Use MPI_Asend() instead of MPI_Send() in CollectData() and copy().
    //
    // Turn off via ParmParse using "fabarray.do_async_sends=0" in inputs file.
    //
    // Default is true.
    //
    static bool do_async_sends;
    //
    // Initialize from ParmParse with "fabarray" prefix.
    //
    static void Initialize ();
    static void Finalize ();
    bool IsInitialized() const;
    void SetInitialized(bool binit);
    //
    // To maximize thread efficiency we now can decompose things like
    // intersections among boxes into smaller tiles. This sets
    // their maximum size.
    //
    static IntVect comm_tile_size;  // communication tile size
    //
    // The number of FabArrays that have been defined during this run
    //
    static int nFabArrays;
    static int NFabArrays() { return nFabArrays; }
    int AllocatedFAPtrID() const  { return aFAPId; }

    struct FPinfo
    {
	FPinfo (const FabArrayBase& srcfa,
		const FabArrayBase& dstfa,
		Box                 dstdomain,
		int                 dstng,
		const BoxConverter& coarsener);
	~FPinfo ();

	long bytes () const;

	BoxArray            ba_crse_patch;
	DistributionMapping dm_crse_patch;
	Array<int>          dst_idxs;
	Array<Box>          dst_boxes;
	//
	BDKey               m_srcbdk;
	BDKey               m_dstbdk;
	Box                 m_dstdomain;
	int                 m_dstng;
	BoxConverter*       m_coarsener;
	//
	int                 m_nuse;
    };

    typedef std::multimap<BDKey,FabArrayBase::FPinfo*> FPinfoCache;
    typedef FPinfoCache::iterator FPinfoCacheIter;

    static FPinfoCache m_TheFillPatchCache;

    static CacheStats m_FPinfo_stats;

    static const FPinfo& TheFPinfo (const FabArrayBase& srcfa,
				    const FabArrayBase& dstfa,
				    Box                 dstdomain,
				    int                 dstng,
				    const BoxConverter& coarsener);

    void flushFPinfo (bool no_assertion=false);

    //
    // parallel copy or add
    //
    enum CpOp { COPY = 0, ADD = 1 };

protected:

    DistributionMapping& ModifyDistributionMap () { return distributionMap; }

    //
    // Returns owenership of fabs. The concept of ownership only applies when UPC++
    // team is used. In that case, each fab is shared by team workers, with one
    // taking the ownership.
    //
    const std::vector<bool>& OwnerShip () const { return ownership; }
    bool isOwner (int li) const { return ownership[li]; }

    //
    // The data ...
    //
    mutable BoxArray    boxarray;
    DistributionMapping distributionMap;
    Array<int>          indexArray;
    std::vector<bool>   ownership;
    int                 n_grow;
    int                 n_comp;
    int                 aFAPId;      // ---- id of currently allocated fab array pointers
    int                 aFAPIdLock;  // ---- lock for resizing sidecars

    mutable BDKey m_bdkey;

    //
    // Tiling
    //
    // We use tile size as the key for the inner map.
    typedef std::map<IntVect, TileArray, IntVect::Compare> TAMap;
    typedef std::map<BDKey, TAMap> TACache;
    //
    static TACache     m_TheTileArrayCache;
    static CacheStats  m_TAC_stats;
    //
    const TileArray* getTileArray (const IntVect& tilesize) const;
    void buildTileArray (const IntVect& tilesize, TileArray& ta) const;
    //
    void flushTileArray (const IntVect& tilesize = IntVect::TheZeroVector(), 
			 bool no_assertion=false) const;
    static void flushTileArrayCache (); // This flushes the entire cache.

    //
    // FillBoundary
    //
    struct FB
    {
        FB (const FabArrayBase& fa, bool cross, const Periodicity& period,
	    bool enforce_periodicity_only);
        ~FB ();

	IndexType           m_typ;
        int                 m_ngrow;
        bool                m_cross;
	bool                m_epo;
	Periodicity         m_period;
        //
        // The cache of local and send/recv per FillBoundary().
        //
	bool                m_threadsafe_loc;
	bool                m_threadsafe_rcv;
        CopyComTagsContainer*      m_LocTags;
        MapOfCopyComTagContainers* m_SndTags;
        MapOfCopyComTagContainers* m_RcvTags;
        std::map<int,int>*         m_SndVols;
        std::map<int,int>*         m_RcvVols;
	//
	int                 m_nuse;
	//
	long bytes () const;
    private:
	void define_fb (const FabArrayBase& fa);
	void define_epo (const FabArrayBase& fa);
    };
    //
    typedef std::multimap<BDKey,FabArrayBase::FB*> FBCache;
    typedef FBCache::iterator FBCacheIter;
    //
    static FBCache    m_TheFBCache;
    static CacheStats m_FBC_stats;
    //
    const FB& getFB (const Periodicity& period, bool cross=false, bool enforce_periodicity_only = false) const;
    //
    void flushFB (bool no_assertion=false) const;       // This flushes its own FB.
    static void flushFBCache (); // This flushes the entire cache.

    //
    // parallel copy or add
    //
    struct CPC
    {
	CPC (const FabArrayBase& dstfa, int dstng,
	     const FabArrayBase& srcfa, int srcng,
	     const Periodicity& period);
	CPC (const BoxArray& dstba, const DistributionMapping& dstdm, 
	     const Array<int>& dstidx, int dstng,
	     const BoxArray& srcba, const DistributionMapping& srcdm, 
	     const Array<int>& srcidx, int srcng,
	     const Periodicity& period, int myproc);
        ~CPC ();

        long bytes () const;	

	BDKey       m_srcbdk;
	BDKey       m_dstbdk;
	int         m_srcng;
	int         m_dstng;
	Periodicity m_period;
	BoxArray    m_srcba;
	BoxArray    m_dstba;
        //
        // The cache of local and send/recv info per FabArray::copy().
        //
	bool        m_threadsafe_loc;
	bool        m_threadsafe_rcv;
        CopyComTagsContainer*      m_LocTags;
        MapOfCopyComTagContainers* m_SndTags;
        MapOfCopyComTagContainers* m_RcvTags;
        std::map<int,int>*         m_SndVols;
        std::map<int,int>*         m_RcvVols;
	//
        int         m_nuse;

    private:
	void define (const BoxArray& ba_dst, const DistributionMapping& dm_dst,
		     const Array<int>& imap_dst,
		     const BoxArray& ba_src, const DistributionMapping& dm_src,
		     const Array<int>& imap_src,
		     int MyProc = ParallelDescriptor::MyProc());
    };
    //
    typedef std::multimap<BDKey,FabArrayBase::CPC*> CPCache;
    typedef CPCache::iterator CPCacheIter;
    //
    static CPCache    m_TheCPCache;
    static CacheStats m_CPC_stats;
    //
    const CPC& getCPC (int dstng, const FabArrayBase& src, int srcng, const Periodicity& period) const;
    // 
    void flushCPC (bool no_assertion=false) const;      // This flushes its own CPC.
    static void flushCPCache (); // This flusheds the entire cache.

    //
    // Keep track of how many FabArrays are built with the same BDKey.
    //
    static std::map<BDKey, int> m_BD_count;
    //
    // clear BD count and caches associated with this BD, if no other is using this BD. 
    // 
    void clearThisBD (bool no_assertion=false);
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
        AllBoxes      = 0x02, // If on, all threads/workers loop over all boxes without tiling.
	                      // This essentially loops over indexMap.
	                      // Note that many functions won't work with this.
        NoTeamBarrier = 0x04, // For Team only. If on, there is no barrier in MFIter dtor.
	SkipInit      = 0x08  // Used by MFGhostIter
    };  // All these flags are off by default.
    //
    // Construct a MFIter.
    //
    // The default is no tiling
    explicit MFIter (const FabArrayBase& fabarray,
		     unsigned char       flags_=0);
    // tiling w/ default size, IntVect FabArrayBase::mfiter_tile_size
    MFIter (const FabArrayBase& fabarray, 
	    bool                do_tiling); 
    // tiling with explicit size and flags
    MFIter (const FabArrayBase& fabarray, 
	    const IntVect&      tilesize,
	    unsigned char       flags_=0);
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
    Box validbox () const;
    //
    // Returns the Box of the FAB at which we currently point.
    //
    Box fabbox () const;
    //
    // Increments iterator to the next tile we own.
    //
    void operator++ () { ++currentIndex; }
    //
    // Is the iterator valid i.e. is it associated with a FAB?
    //
    bool isValid () { return currentIndex < endIndex; }
    //
    // The index into the underlying BoxArray of the current FAB.
    //
    int index () const { return (*index_map)[currentIndex]; }
    //
    // local index into the vector of fab pointers, m_fabs_v
    // When AllBoxes is on, local_index_map is a nullptr and local index is current index.
    //
    int LocalIndex () const { return local_index_map ? (*local_index_map)[currentIndex] : currentIndex; }
    //
    // Constant reference to FabArray over which we're iterating.
    //
    const FabArrayBase& theFabArrayBase () const { return fabArray; }

protected:

    const FabArrayBase& fabArray;

    IntVect tile_size;

    unsigned char flags;
    int           currentIndex;
    int           beginIndex;
    int           endIndex;
    IndexType     typ;

    const Array<int>* index_map;
    const Array<int>* local_index_map;
    const Array<Box>* tile_array;

    void Initialize ();
};

/*
 * Iterate over ghost cells.  Lots of MFIter functons do not work.
 */
class MFGhostIter
    :
    public MFIter
{
public:
    explicit MFGhostIter (const FabArrayBase& fabarray);
private:
    void Initialize ();
    FabArrayBase::TileArray lta;
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

    FabArray (const BoxArray& bxs,
              int             nvar,
              int             ngrow,
	      ParallelDescriptor::Color color);
    //
    // The destructor -- deletes all FABs in the array.
    //
    virtual ~FabArray ();
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
		 const IntVect&  nodal = IntVect::TheZeroVector(),
		 ParallelDescriptor::Color color = ParallelDescriptor::DefaultColor());

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
	       const Periodicity&   period = Periodicity::NonPeriodic(),
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
	       const Periodicity&   period = Periodicity::NonPeriodic(),
               CpOp                 op = FabArrayBase::COPY);
    void copy (const FabArray<FAB>& src,
               int                  src_comp,
               int                  dest_comp,
               int                  num_comp,
	       int                  src_nghost,
	       int                  dst_nghost,
	       const Periodicity&   period = Periodicity::NonPeriodic(),
               CpOp                 op = FabArrayBase::COPY);

    //
    // In the following copyTo functions, the destination FAB is identical on each process!!
    //
    // Copies the values contained in the intersection of the
    // valid + nghost region of this FabArray with the FAB dest into dest.
    //
    void copyTo (FAB& dest,
		 int  nghost = 0) const;
    //
    // Copies the values contained in the intersection of the
    // valid + nghost region of this FabArray with the FAB dest and the Box
    // subbox into that subregion of dest.
    //
    void copyTo (FAB&       dest,
                const Box& subbox,
		 int        nghost = 0) const;
    //
    // Copies the values contained in the intersection of the
    // num_comp component valid + nghost region of this FabArray, starting at
    // component src_comp, with the FAB dest into dest, starting at
    // component dest_comp in dest.
    //
    void copyTo (FAB& dest,
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
    void copyTo (FAB&       dest,
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
    // the FabArray.  If cross=true, corner cells are not filled.
    // If the length of periodic is provided, periodic boundaries are
    // also filled.  Note that FabArray itself does not contains
    // any periodicity information.
    // FillBoundary expects that its cell-centered version of its BoxArray 
    // is non-overlapping.
    //
    void FillBoundary (bool cross = false);
    void FillBoundary (const Periodicity& period, bool cross = false);
    //
    // Same as FillBoundary(), but only copies ncomp components starting at scomp.
    //
    void FillBoundary (int scomp, int ncomp, bool cross = false);
    void FillBoundary (int scomp, int ncomp, const Periodicity& period, bool cross = false);

    void FillBoundary_nowait (bool cross = false);
    void FillBoundary_nowait (const Periodicity& period, bool cross = false);
    void FillBoundary_nowait (int scomp, int ncomp, bool cross = false);
    void FillBoundary_nowait (int scomp, int ncomp, const Periodicity& period, bool cross = false);
    void FillBoundary_finish ();

    // Fill cells outside periodic domains with their corresponding cells inside
    // the domain.  Ghost cells are treated the same as valid cells.  The BoxArray
    // is allowed to be overlapping.
    void EnforcePeriodicity (const Periodicity& period);
    void EnforcePeriodicity (int scomp, int ncomp, const Periodicity& period);

    // covered  : ghost cells covered by valid cells of this FabArray
    //            (including periodically shifted valid cells)
    // uncovered: ghost cells not covered by valid cells
    //            (including ghost cells outside periodic boundaries)
    // physbnd  : boundary cells outside the domain (excluding periodic boundaries)
    // interior : interior cells (i.e., valid cells)
    void BuildMask (const Box& phys_domain, const Periodicity& period,
		    value_type covered, value_type uncovered, 
		    value_type physbnd, value_type interior);

    //
    // Move FABs in this FabArray to different MPI ranks.
    //

    struct FABMoves {
      int distMapIndex, fromRank, toRank, seqNum;
    };

    int MoveFabs (const Array<int> &newDistMapArray);
    static void MoveAllFabs (const Array<int> &newDistMapArray);
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
    //
    // Helper function for define().
    //
    void defineDoit (const BoxArray&            bxs,
                     int                        nvar,
                     int                        ngrow,
                     FabAlloc                   mem_mode,
                     const DistributionMapping* dm,
		     const IntVect&             nodal,
		     ParallelDescriptor::Color  color);
    //
    // The data.
    //
    std::vector<FAB*> m_fabs_v;

    // ---- currently active fabarrays:  <ngrids to find distmap, <aFAPId, FabArray *> >
    static std::map<int, std::map<int, FabArray<FAB> *> > allocatedFAPointers;

    // ---- set of ids of fabarrays created with Fab_noalloc
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
	    if (alloc)
		BoxLib::update_fab_stats(-n_points, -n_values, sizeof(value_type));
	}
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
    //
    // These are disallowed.
    //
    FabArray (const FabArray<FAB>&);
    FabArray<FAB>& operator= (const FabArray<FAB>&);
    //
    // This is used locally in all define functions.
    //
    void AllocFabs ();

    void FBEP_nowait (int scomp, int ncomp, const Periodicity& period, bool cross,
		      bool enforce_periodicity_only = false);

public:
    // Data used in non-blocking FillBoundary
    bool fb_cross, fb_epo;
    int fb_scomp, fb_ncomp;
    Periodicity fb_period;

    //
    value_type*        fb_the_recv_data;
    Array<int>         fb_recv_from;
    Array<value_type*> fb_recv_data;
    Array<MPI_Request> fb_recv_reqs;
#ifdef BL_USE_MPI3
    Array<MPI_Aint>    fb_recv_disp;
#endif
    //
    Array<value_type*> fb_send_data;
    Array<MPI_Request> fb_send_reqs;

    // Data used in non-blocking FillPeriodicBoundary
    bool fpb_corners;
    int fpb_scomp, fpb_ncomp;
    //
    value_type*        fpb_the_recv_data;
    Array<int>         fpb_recv_from;
    Array<value_type*> fpb_recv_data;
    Array<MPI_Request> fpb_recv_reqs;
#ifdef BL_USE_MPI3
    Array<MPI_Aint>    fpb_recv_disp;
#endif
    //
    Array<value_type*> fpb_send_data;
    Array<MPI_Request> fpb_send_reqs;
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
#if BL_USE_UPCXX
    static struct FabArrayCopyDescriptor<FAB>* CurrentCollectDataObj;
    static std::vector<void*> CollectData_send_data_buffers;

    static void CollectData_AM_handler(uint32_t src_rank,
                                       void *payload,
                                       size_t nbytes,
                                       void *target_buffer,
                                       void *target_event,
                                       int numboxes,
                                       int not_used);
#endif

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

    int CurrentNFabArrays () const { return fabArrays.size(); }

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
FabArrayBase::PostRcvs (const std::map<int,int>&               m_RcvVols,
                        T*&                                    the_recv_data,
                        Array<T*>&                             recv_data,
                        Array<int>&                            recv_from,
                        Array<MPI_Request>&                    recv_reqs,
                        int                                    ncomp,
                        int                                    SeqNum,
			MPI_Comm                               comm)
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
        recv_reqs.push_back(ParallelDescriptor::Arecv(recv_data.back(),N,it->first,SeqNum,comm).req());

        Offset += N;
    }
}

#ifdef BL_USE_MPI3
template<typename T>
void
FabArrayBase::PostRcvs_MPI_Onesided (
			const std::map<int,int>&               m_RcvVols,
                        T*&                                    the_recv_data,
                        Array<T*>&                             recv_data,
                        Array<int>&                            recv_from,
                        Array<MPI_Request>&                    recv_reqs,
                        Array<MPI_Aint> &                      recv_disp,
                        int                                    ncomp,
			int                                    SeqNum,
			MPI_Win &                              win)
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

    MPI_Win_attach(win, the_recv_data, TotalRcvsVolume*sizeof(T));

    int Offset = 0;
    int i = 0;

    recv_disp.resize(m_RcvVols.size());  /// Must resize it now, otherwise push_back may reallocate the vector

    for (std::map<int,int>::const_iterator it = m_RcvVols.begin(),
             End = m_RcvVols.end();
         it != End;
         ++it)
    {
        const int N = it->second*ncomp;

        BL_ASSERT(N < std::numeric_limits<int>::max());

        recv_data.push_back(&the_recv_data[Offset]);
        recv_from.push_back(it->first);
	MPI_Get_address(&the_recv_data[Offset], &recv_disp[i]);
        recv_reqs.push_back(ParallelDescriptor::Asend(&recv_disp[i],1,it->first,SeqNum).req());

	i++;
        Offset += N;
    }
}
#endif   // BL_USE_MPI3

template<typename T>
void
FabArrayBase::WaitForAsyncSends (int                 N_snds,
                                 Array<MPI_Request>& send_reqs,
                                 Array<T*>&          send_data,
                                 Array<MPI_Status>&  stats)
{
#ifdef BL_USE_MPI
    BL_ASSERT(N_snds > 0);

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
    BL_ASSERT(mfi.LocalIndex() < indexArray.size());
    return *m_fabs_v[mfi.LocalIndex()];
}

template <class FAB>
FAB&
FabArray<FAB>::operator[] (const MFIter& mfi)
{
    BL_ASSERT(mfi.LocalIndex() < indexArray.size());
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
    clearThisBD();

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
    : shmem()
{
    m_FA_stats.recordBuild();
}

template <class FAB>
FabArray<FAB>::FabArray (const BoxArray& bxs,
                         int             nvar,
                         int             ngrow,
                         FabAlloc        alloc,
			 const IntVect&  nodal)
    : shmem()
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
    : shmem()
{
    m_FA_stats.recordBuild();
    define(bxs,nvar,ngrow,dm,alloc,nodal);
}

template <class FAB>
FabArray<FAB>::FabArray (const BoxArray& bxs,
                         int             nvar,
                         int             ngrow,
			 ParallelDescriptor::Color color)
    : shmem()
{
    m_FA_stats.recordBuild();
    define(bxs,nvar,ngrow,Fab_allocate,IntVect::TheZeroVector(),color);
}

template <class FAB>
FabArray<FAB>::~FabArray ()
{
    m_FA_stats.recordDelete();

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
FabArray<FAB>::defineDoit (const BoxArray&            bxs,
                           int                        nvar,
                           int                        ngrow,
                           FabAlloc                   alloc,
                           const DistributionMapping* dm,
			   const IntVect&             nodal,
			   ParallelDescriptor::Color  color)
{
    BL_ASSERT(ngrow >= 0);
    BL_ASSERT(boxarray.size() == 0);
    indexArray.clear();
    ownership.clear();
    n_grow = ngrow;
    n_comp = nvar;

    boxarray = bxs;

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
        distributionMap.define(boxarray,ParallelDescriptor::NProcs(),color);
    }
    else
    {
        BL_ASSERT(dm->ProcessorMap().size() == bxs.size()+1);

        distributionMap = *dm;
    }

    int myProc = ParallelDescriptor::MyProc();

    for(int i = 0, N = boxarray.size(); i < N; ++i) {
	if (ParallelDescriptor::sameTeam(distributionMap[i])) {
	// If Team is not used (i.e., team size == 1), distributionMap[i] == myProc
	    indexArray.push_back(i);
	    ownership.push_back(myProc == distributionMap[i]);
	}
    }

    addThisBD();

    if(alloc == Fab_allocate) {
      AllocFabs();
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

    if(alloc == Fab_allocate) {
      noallocFAPIds.erase(aFAPId);
    } else {
      noallocFAPIds.insert(aFAPId);
    }

#ifdef BL_USE_TEAM
    if (alloc == Fab_allocate) ParallelDescriptor::MyTeam().MemoryBarrier();
#endif
}

template <class FAB>
void
FabArray<FAB>::define (const BoxArray& bxs,
                       int             nvar,
                       int             ngrow,
                       FabAlloc        alloc,
		       const IntVect&  nodal,
		       ParallelDescriptor::Color color)
{
    defineDoit(bxs,nvar,ngrow,alloc,0,nodal,color);
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
    defineDoit(bxs,nvar,ngrow,alloc,&dm,nodal,ParallelDescriptor::DefaultColor());
}

template <class FAB>
void
FabArray<FAB>::AllocFabs ()
{
    const int n = indexArray.size();
    const int nworkers = ParallelDescriptor::TeamSize();
    shmem.alloc = (nworkers > 1);

    m_fabs_v.reserve(n);

    for (int i = 0; i < n; ++i)
    {
	int K = indexArray[i];
        const Box& tmp = fabbox(K);
	bool alloc = !shmem.alloc;
	m_fabs_v.push_back(new FAB(tmp, n_comp, alloc, shmem.alloc));
    }
    
#ifdef BL_USE_TEAM
    if (shmem.alloc)
    {
	const int teamlead = ParallelDescriptor::MyTeamLead();

	shmem.n_values = 0;
	shmem.n_points = 0;
	Array<long> offset(n,0);
	Array<long> nextoffset(nworkers,-1);
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
	Array<value_type*> dps;

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

	BoxLib::Abort("BaseFab::define: to allocate shared memory, either USE_UPCXX or USE_MPI3 must be true");

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

	BoxLib::update_fab_stats(shmem.n_points, shmem.n_values, sizeof(value_type));
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
    if (n_comp == 0)
        n_comp = elem->nComp();

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
                       FAB* elem)
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

    if (m_fabs_v.size() == 0) {
	m_fabs_v.resize(indexArray.size());
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
		     const Periodicity&   period,
                     CpOp                 op)
{
    BL_PROFILE("FabArray::copy()");

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

    const CPC& thecpc = getCPC(dnghost, src, snghost, period);

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
                get(tag.dstIndex).copy(src[tag.srcIndex],tag.sbox,scomp,tag.dbox,dcomp,ncomp);
            }
            else
            {
                get(tag.dstIndex).plus(src[tag.srcIndex],tag.sbox,tag.dbox,scomp,dcomp,ncomp);
            }
        }

        return;
    }

#ifdef BL_USE_MPI

#if defined(BL_USE_UPCXX)
    ParallelDescriptor::Mode.set_upcxx_mode();
#endif

#if !defined(BL_USE_MPI3)
    BL_ASSERT(!ParallelDescriptor::MPIOneSided());
#endif

    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    int SeqNum;
    {
	ParallelDescriptor::Color src_color = src.color();
	ParallelDescriptor::Color dst_color = this->color();
	if (src_color == ParallelDescriptor::DefaultColor() ||
	    dst_color == ParallelDescriptor::DefaultColor() ||
	    src_color != dst_color) {
	    // If any of the two FabArrays is in the global communicator, 
	    // or if the two have different color,
	    // all processes are here.
	    SeqNum = ParallelDescriptor::SeqNum();
	} else { // The two have the same non-default color.
	    if (ParallelDescriptor::SubCommColor() == src_color) {
		SeqNum = ParallelDescriptor::SubSeqNum();
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
    if (ParallelDescriptor::MPIOneSided())
	MPI_Comm_group(ParallelDescriptor::Communicator(this->color()), &tgroup);
#endif

    //
    // Send/Recv at most MaxComp components at a time to cut down memory usage.
    //
    int NCompLeft = ncomp;

    for (int ipass = 0, SC = scomp, DC = dcomp; ipass < ncomp; )
    {
        const int NC = std::min(NCompLeft,FabArrayBase::MaxComp);

        Array<int>         recv_from;
        Array<value_type*> recv_data;
        Array<MPI_Request> recv_reqs;
#ifdef BL_USE_MPI3
	Array<MPI_Aint>    recv_disp;
#endif
        //
        // Post rcvs. Allocate one chunk of space to hold'm all.
        //
        value_type* the_recv_data = 0;

	if (N_rcvs > 0) {
#ifdef BL_USE_UPCXX
	    FabArrayBase::PostRcvs_PGAS(*thecpc.m_RcvVols,the_recv_data,
					recv_data,recv_from,NC,SeqNum,&BLPgas::cp_recv_event);
#else
	    if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
		FabArrayBase::PostRcvs_MPI_Onesided(*thecpc.m_RcvVols,the_recv_data,
						    recv_data,recv_from,recv_reqs,
						    recv_disp,NC,SeqNum,ParallelDescriptor::cp_win);
		MPI_Group_incl(tgroup, recv_from.size(), recv_from.dataPtr(), &rgroup);
		MPI_Win_post(rgroup, 0, ParallelDescriptor::cp_win);
#endif
	    } else {
		FabArrayBase::PostRcvs(*thecpc.m_RcvVols,the_recv_data,
				       recv_data,recv_from,recv_reqs,NC,SeqNum);
	    }
#endif
	}

	//
	// Post send's
	// 
	Array<value_type*>                 send_data;
	Array<int>                         send_N;
	Array<int>                         send_rank;
	Array<MPI_Request>                 send_reqs;
	Array<const CopyComTagsContainer*> send_cctc;

	if (N_snds > 0)
	{
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

		value_type* data = static_cast<value_type*>
#ifdef BL_USE_UPCXX
		    (BLPgas::alloc(N*sizeof(value_type)));
#else
		    (BoxLib::The_Arena()->alloc(N*sizeof(value_type)));
#endif
 
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
		    const Box& bx = it->sbox;
		    src[it->srcIndex].copyToMem(bx,SC,NC,dptr);
		    const int Cnt = bx.numPts()*NC;
		    dptr += Cnt;
		}
	    }

#ifdef BL_USE_UPCXX
	    
	    BLPgas::cp_send_counter = 0;

	    for (int i=0; i<N_snds; ++i) {
		BLPgas::Send(upcxx::global_ptr<void>((void *)send_data[i], upcxx::myrank()),
			     send_rank[i],
			     send_N[i]*sizeof(value_type),
			     SeqNum,
			     &BLPgas::cp_send_event,
			     &BLPgas::cp_send_counter);
	    }
	    
	    // Need to make sure at least half of the sends have been started
	    while (BLPgas::cp_send_counter < N_snds)
		upcxx::advance();

#else
	    
	    if (ParallelDescriptor::MPIOneSided())
	    {
#if defined(BL_USE_MPI3)
		Array<MPI_Aint> send_disp;
		send_reqs.reserve(N_snds);
		send_disp.resize (N_snds);
		
		for (int i=0; i<N_snds; ++i) {
		    send_reqs.push_back(ParallelDescriptor::Arecv
					(&send_disp[i],1,send_rank[i],SeqNum).req());
		}
		
		MPI_Group_incl(tgroup, send_rank.size(), send_rank.dataPtr(), &sgroup);
		MPI_Win_start(sgroup,0,ParallelDescriptor::cp_win);
		
		int send_counter = 0;	
		while (send_counter < N_snds) {
		    MPI_Status status;
		    int index;
		    
		    MPI_Waitany(N_snds, send_reqs.dataPtr(), &index, &status);
		    
		    BL_ASSERT(status.MPI_TAG == SeqNum);
		    BL_ASSERT(status.MPI_SOURCE == send_rank[index]);
		    
		    MPI_Put(send_data[index], send_N[index]*sizeof(value_type), MPI_CHAR, send_rank[index],
			    send_disp[index], send_N[index]*sizeof(value_type), MPI_CHAR, 
			    ParallelDescriptor::cp_win);
		    
		    send_counter++;
		}
#endif
	    } else {
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
#pragma omp parallel
#endif
	    ParallelDescriptor::team_for(0, N_locs, [&] (int j) 
            {
		const CopyComTag& tag = (*thecpc.m_LocTags)[j];
		
		if (op == FabArrayBase::COPY)
		{
		    get(tag.dstIndex).copy(src[tag.srcIndex],tag.sbox,SC,tag.dbox,DC,NC);
		}
		else
		{
		    get(tag.dstIndex).plus(src[tag.srcIndex],tag.sbox,tag.dbox,SC,DC,NC);
		}
	    });
#endif	    
	}
	else 
	{
#ifdef _OPENMP
#pragma omp parallel for if (thecpc.m_threadsafe_loc)
#endif
	    for (int j=0; j<N_locs; ++j)
	    {
		const CopyComTag& tag = (*thecpc.m_LocTags)[j];
		
		if (op == FabArrayBase::COPY)
		{
		    get(tag.dstIndex).copy(src[tag.srcIndex],tag.sbox,SC,tag.dbox,DC,NC);
		}
		else
		{
		    get(tag.dstIndex).plus(src[tag.srcIndex],tag.sbox,tag.dbox,SC,DC,NC);
		}
	    }
	}

	//
	//  wait and unpack
	//
#ifdef BL_USE_UPCXX
        if (N_rcvs > 0) BLPgas::cp_recv_event.wait();
#else
	if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
	    if (N_snds > 0) MPI_Win_complete(ParallelDescriptor::cp_win);
	    if (N_rcvs > 0) MPI_Win_wait    (ParallelDescriptor::cp_win);
#endif
	} else {
	    if (N_rcvs > 0) {
		Array<MPI_Status> stats(N_rcvs);
		BL_MPI_REQUIRE( MPI_Waitall(N_rcvs, recv_reqs.dataPtr(), stats.dataPtr()) );
	    }
	}
#endif	

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
			const Box& bx  = it->dbox;
			const int  Cnt = bx.numPts()*NC;
			
			if (op == FabArrayBase::COPY)
			{
			    get(it->dstIndex).copyFromMem(bx,DC,NC,dptr);
			}
			else
			{
			    fab.resize(bx,NC);
			    memcpy(fab.dataPtr(), dptr, Cnt*sizeof(value_type));
			    get(it->dstIndex).plus(fab,bx,bx,0,DC,NC);
			}
			
			dptr += Cnt;
		    }
		}
	    }

#ifdef BL_USE_UPCXX
	    BLPgas::free(the_recv_data);
#else
	    if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
		MPI_Win_detach(ParallelDescriptor::cp_win, the_recv_data);
		BoxLib::The_Arena()->free(the_recv_data);
		recv_disp.clear();
#endif
	    } else {
		BoxLib::The_Arena()->free(the_recv_data);
	    }
#endif
	    recv_from.clear();
	    recv_data.clear();
	    recv_reqs.clear();
	}
	
        if (N_snds > 0) {
#ifdef  BL_USE_UPCXX
	    FabArrayBase::WaitForAsyncSends_PGAS(N_snds,send_data,
					         &BLPgas::cp_send_event,
					         &BLPgas::cp_send_counter);
#else
	    if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
		for (int i = 0; i < N_snds; ++i)
		    BoxLib::The_Arena()->free(send_data[i]);
#endif
	    } else {
		if (FabArrayBase::do_async_sends && ! thecpc.m_SndTags->empty()) {
		    Array<MPI_Status> stats;
		    FabArrayBase::WaitForAsyncSends(thecpc.m_SndTags->size(),send_reqs,send_data,stats);
		}
	    }
#endif
            send_data.clear();
	    send_reqs.clear();
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
FabArray<FAB>::copy (const FabArray<FAB>& src,
                     int                  scomp,
                     int                  dcomp,
                     int                  ncomp,
		     const Periodicity&   period,
                     CpOp                 op)
{
    copy(src,scomp,dcomp,ncomp,0,0,period,op);
}

template <class FAB>
void
FabArray<FAB>::copy (const FabArray<FAB>& src, const Periodicity& period, CpOp op)
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
	    const Box& bx = BoxLib::grow(boxarray[j],nghost);
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

    bool putInCache(false);
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
    Array<int> srcPMap, destPMap, srcPMapAll, destPMapAll;
    DistributionMapping srcDM, destDM, srcDMAll, destDMAll;
    Array<int> srcIndexArray, destIndexArray;
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
      srcPMapAll[srcPMapAll.size()-1] = myProcBoth;  // ---- set the sentinel
    }
    if(isDest && myProcBoth == bcastRootDest) {
      destBoxArray = dest->boxArray();
      destDM = dest->DistributionMap();
      destPMap = destDM.ProcessorMap();
      destPMapAll = DistributionMapping::TranslateProcMap(destPMap, groupBoth, groupDest);
      destPMapAll[destPMapAll.size()-1] = myProcBoth;  // ---- set the sentinel
    }

    BoxLib::BroadcastBoxArray(srcBoxArray,  myProcBoth, bcastRootSrc,  commBoth);

    BoxLib::BroadcastBoxArray(destBoxArray, myProcBoth, bcastRootDest, commBoth);
    BL_ASSERT(destBoxArray[0].ixType() == srcBoxArray[0].ixType());


    BoxLib::BroadcastArray(srcPMapAll,  myProcBoth, bcastRootSrc,  commBoth);
    BoxLib::BroadcastArray(destPMapAll, myProcBoth, bcastRootDest, commBoth);

    srcDMAll.define(srcPMapAll, putInCache);
    destDMAll.define(destPMapAll, putInCache);

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

        FabArrayBase::PostRcvs(*thecpc.m_RcvVols, the_recv_data, recv_data,
	                       recv_from, recv_reqs, NC, SeqNum, commBoth);

	//
	// Post sends
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
                const Box& bx = it->sbox;
                (*src)[it->srcIndex].copyToMem(bx,SC,NC,dptr);
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
				    (send_data[j],send_N[j],send_rank[j],SeqNum, commBoth).req());
            }
	} else {
	    for (int j=0; j<N_snds; ++j)
	    {
                ParallelDescriptor::Send(send_data[j],send_N[j],send_rank[j],SeqNum, commBoth);
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
		    const Box& bx  = it->dbox;
		    const int  Cnt = bx.numPts()*NC;
		    
		    if (op == FabArrayBase::COPY)
		    {
			dest->get(it->dstIndex).copyFromMem(bx,DC,NC,dptr);
		    }
		    else
		    {
			fab.resize(bx,NC);
			memcpy(fab.dataPtr(), dptr, Cnt*sizeof(value_type));
			
			dest->get(it->dstIndex).plus(fab,bx,bx,0,DC,NC);
		    }
		    
		    dptr += Cnt;
		}
	    }
	}
	}
	
        BoxLib::The_Arena()->free(the_recv_data);
	
        if (FabArrayBase::do_async_sends && ! thecpc.m_SndTags->empty()) {
            FabArrayBase::WaitForAsyncSends(thecpc.m_SndTags->size(),send_reqs,send_data,stats);
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
void
FabArray<FAB>::setVal (value_type val)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter fai(*this,true); fai.isValid(); ++fai)
    {
	const Box& bx = fai.growntilebox();
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
	const Box& bx = fai.growntilebox(nghost);
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
      BoxLib::USleep(ParallelDescriptor::MyProcAll());
      std::cout << ParallelDescriptor::MyProcAll() << ":: **** ERROR:  allocFAPSize != "
                << "FabArray<FAB>::allocatedFAPointers.size() : " << allocFAPSize
	        << "  " << FabArray<FAB>::allocatedFAPointers.size() << std::endl;
      BoxLib::Abort("**** Error in FabArray<FAB>::CheckFAPointers():  bad allocFAPSize.");
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
        BoxLib::USleep(ParallelDescriptor::MyProcAll());
        std::cout << ParallelDescriptor::MyProcAll() << ":: **** ERROR:  fapCMapSize != "
	          << "faPtrCachedMap.size() : " << fapCMapSize << "  "
		  << faPtrCachedMap.size() << std::endl;
        BoxLib::Abort("**** Error in FabArray<FAB>::CheckFAPointers():  bad fapCMapSize.");
      }
    }

    int distmapNGrids(afapIter->first-1);
    ParallelDescriptor::Bcast(&distmapNGrids, 1, iopNumber);
    if(distmapNGrids != afapIter->first-1) {
      if(abortOnError) {
        BoxLib::USleep(ParallelDescriptor::MyProcAll());
        std::cout << ParallelDescriptor::MyProcAll() << ":: **** ERROR:  distmapNGrids != "
	          << "afapIter->first-1 : " << distmapNGrids << "  " << afapIter->first-1 << std::endl;
        BoxLib::Abort("**** Error in FabArray<FAB>::CheckFAPointers():  bad distmapNGrids.");
      }
    }

    Array<int> fapId(faPtrCachedMap.size()), dmId(faPtrCachedMap.size());
    Array<int> nDM(faPtrCachedMap.size()), fabsAlloced(faPtrCachedMap.size());
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
    Array<int> fapIdCheck, dmIdCheck, nDMCheck, fabsAllocedCheck;
    if(ParallelDescriptor::IOProcessor()) {
      fapIdCheck = fapId;
      dmIdCheck = dmId;
      nDMCheck = nDM;
      fabsAllocedCheck = fabsAlloced;
    }
    BoxLib::BroadcastArray(fapIdCheck, ParallelDescriptor::MyProc(), iopNumber,
                           ParallelDescriptor::Communicator());
    BoxLib::BroadcastArray(dmIdCheck, ParallelDescriptor::MyProc(), iopNumber,
                           ParallelDescriptor::Communicator());
    BoxLib::BroadcastArray(nDMCheck, ParallelDescriptor::MyProc(), iopNumber,
                           ParallelDescriptor::Communicator());
    BoxLib::BroadcastArray(fabsAllocedCheck, ParallelDescriptor::MyProc(), iopNumber,
                           ParallelDescriptor::Communicator());
    for(int i(0); i < fapId.size(); ++i) {
      if(fapId[i] != fapIdCheck[i]) {
        if(abortOnError) {
          std::cout << ParallelDescriptor::MyProcAll() << ":: **** ERROR:  fapId[" << i
	            << "] != fapIdCheck[" << i << "] : "
                    << fapId[i] << "  " << fapIdCheck[i] << std::endl;
          BoxLib::Abort("**** Error in FabArray<FAB>::CheckFAPointers():  bad fapId.");
        }
      }
    }
    for(int i(0); i < dmId.size(); ++i) {
      if(dmId[i] != dmIdCheck[i]) {
        if(abortOnError) {
          std::cout << ParallelDescriptor::MyProcAll() << ":: **** ERROR:  dmId[" << i
	            << "] != dmIdCheck[" << i << "] : "
                    << dmId[i] << "  " << dmIdCheck[i] << std::endl;
          BoxLib::Abort("**** Error in FabArray<FAB>::CheckFAPointers():  bad dmId.");
        }
      }
    }
    for(int i(0); i < nDM.size(); ++i) {
      if(nDM[i] != nDMCheck[i]) {
        if(abortOnError) {
          std::cout << ParallelDescriptor::MyProcAll() << ":: **** ERROR:  nDM[" << i
	            << "] != nDMCheck[" << i << "] : "
                    << nDM[i] << "  " << nDMCheck[i] << std::endl;
          BoxLib::Abort("**** Error in FabArray<FAB>::CheckFAPointers():  bad nDM.");
        }
      }
    }
    for(int i(0); i < fabsAlloced.size(); ++i) {
      if(fabsAlloced[i] != fabsAllocedCheck[i]) {
        if(abortOnError) {
          std::cout << ParallelDescriptor::MyProcAll() << ":: **** ERROR:  fabsAlloced[" << i
	            << "] != fabsAllocedCheck[" << i << "] : "
                    << fabsAlloced[i] << "  " << fabsAllocedCheck[i] << std::endl;
          BoxLib::Abort("**** Error in FabArray<FAB>::CheckFAPointers():  bad fabsAlloced.");
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
      BoxLib::Abort("**** Error in FabArray<FAB>::CheckFAPointers():  bad pLocked.");
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
  DistributionMapping::CacheStats(std::cout, ParallelDescriptor::MyProcAll());
  std::cout << "|||| |||| |||| _out PrintFAPointers: ================ " << std::endl;
}


template <class FAB>
void
FabArray<FAB>::AddProcsToComp (int ioProcNumSCS, int ioProcNumAll,
                               int scsMyId, MPI_Comm scsComm)
{
  flushTileArrayCache();

  // ---- BoxArrays
  BoxLib::BroadcastBoxArray(boxarray, scsMyId, ioProcNumSCS, scsComm);

  // ---- DistributionMapping
  int sentinelProc(ParallelDescriptor::MyProcComp());
  BoxLib::BroadcastDistributionMapping(distributionMap, sentinelProc,
                                       scsMyId, ioProcNumSCS, scsComm, true);

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
  Array<int> allIntVects;
  if(scsMyId == ioProcNumSCS) {
    for(int i(0); i < BL_SPACEDIM; ++i)    { allIntVects.push_back(mfiter_tile_size[i]); }
    for(int i(0); i < BL_SPACEDIM; ++i)    { allIntVects.push_back(comm_tile_size[i]); }
  }
  BoxLib::BroadcastArray(allIntVects, scsMyId, ioProcNumSCS, scsComm);
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
  Array<int> setArray;
  if(scsMyId == ioProcNumSCS) {
    std::set<int>::iterator it;
    for(it = noallocFAPIds.begin(); it != noallocFAPIds.end(); ++it) {
      setArray.push_back(*it);
    }
  }
  BoxLib::BroadcastArray(setArray, scsMyId, ioProcNumSCS, scsComm);
  if(scsMyId != ioProcNumSCS) {
    for(int i(0); i < setArray.size(); ++i) {
      noallocFAPIds.insert(setArray[i]);
    }
  }

  // ---- m_bdkey
  if(scsMyId != ioProcNumSCS) {
    m_bdkey =  getBDKey();
  }

  // ---- static TACache     m_TheTileArrayCache;     // ---- cleared with flushTileArrayCache();
  // ---- static CacheStats  m_TAC_stats;             // ---- leave as is for now
  // ---- m_BD_count                                  // ---- leave as is for now
  // ---- m_FA_stats                                  // ---- leave as is for now


  // ---- unlock
  aFAPIdLock = 0;
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
      FabArray<FAB> *faPtr = it->second;
      if(faPtr->ok() == false) {
        BoxLib::Abort("**** Error 0 in MoveAllFabs:  faPtr not ok");
      }

      faPtr->flushFPinfo();
      nFabsMoved = faPtr->MoveFabs(newDistMapArray);  // ---- just use the last return value

      lastFAPtr = faPtr;

    }

    ParallelDescriptor::ReduceIntSum(nFabsMoved);

    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "FabArray<FAB>::MoveAllFabs:  moved " << nFabsMoved
	        << " FAB(s) for each FabArray." << std::endl;
    }

    // ---- flush caches
    DistributionMapping::FlushCache();
    flushCPCache();
    flushFBCache();
    flushTileArrayCache();

    if(lastFAPtr != 0) {
      lastFAPtr->ModifyDistributionMap().ReplaceCachedProcessorMap(newDistMapArray);
    }


    // ---- check the moved fabarray
    for(typename std::map<int, FabArray<FAB> *>::iterator it = faPtrCachedMap.begin();
        it != faPtrCachedMap.end(); ++it)
    {
      FabArray<FAB> *faPtr = it->second;
      if(faPtr->ok() == false) {
        BoxLib::Abort("**** Error 2 in MoveAllFabs:  faPtr not ok");
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
      std::cout << "**** idm ndm[idm] np = " << newDistMapArray[idm] << "  "
                << ParallelDescriptor::NProcs() << std::endl;
      BoxLib::Abort("**** Error:  bad newDistMap:2");
    }
  }
  if(newDistMapArray == distributionMap.ProcessorMap()) {
    return 0;
  }

  // ---- try to resolve fabarrays that were created with Fab_noallocate
  std::set<int>::iterator it = noallocFAPIds.find(aFAPId);
  bool fabsNotAllocated(it != noallocFAPIds.end());
  if(fabsNotAllocated) {
    if(ok()) {  // ---- this fabarray was created with Fab_noallocate then allocated later
      noallocFAPIds.erase(aFAPId);
    } else {
      return 0;
    }
  } else {
    BL_ASSERT(ok());
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
    if(indexArray.size() != m_fabs_v.size()) {
      BoxLib::Abort("**** Error:  indexArray size != m_fabs_v.size()");
    }
    for(int iim(0); iim < indexArray.size(); ++iim) {
      tIFiter = tempIndexFABs.find(indexArray[iim]);
      if(tIFiter == tempIndexFABs.end()) {
	tempIndexFABs.insert(std::pair<int, FAB *>(indexArray[iim], m_fabs_v[iim]));
      } else {
        BoxLib::Abort("**** Error:  index not in indexArray.");
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
  indexArray.clear();
  ownership.clear();
  m_fabs_v.clear();
  for(tIFiter = tempIndexFABs.begin(); tIFiter != tempIndexFABs.end(); ++tIFiter) {
    indexArray.push_back(tIFiter->first);
    ownership.push_back(myProc == newDistMapArray[tIFiter->first]);
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

    {  // ---- check distmap
      const DistributionMapping &dm = fabarray->DistributionMap();
      if(fabarray != 0) {
        dm.Check();
      }
    }

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
    const int myProc = ParallelDescriptor::MyProc();

    FabArray<FAB>* fabArray = fabArrays[fabarrayid.Id()];

    fabArray->DistributionMap().Check();

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
	if(remoteProc >= ParallelDescriptor::NProcs()) {
	  std::cout << ParallelDescriptor::MyProcAll() << ":: _in AddBoxDoIt:  nProcs remoteProc = "
	            << ParallelDescriptor::NProcs() << "  " << remoteProc << std::endl;
	  BoxLib::Abort("Bad remoteProc");
	}
        fcd->fillBoxId     = nextFillBoxId;
        fcd->subBox        = intersect;
        fcd->myProc        = myProc;
        fcd->copyFromProc  = remoteProc;
        fcd->copyFromIndex = faindex;
        fcd->srcComp       = srccomp;
        fcd->destComp      = destcomp;
        fcd->nComp         = numcomp;

        if (ParallelDescriptor::sameTeam(remoteProc))
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
            fabComTag.procThatNeedsData = myProc;
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

#if BL_USE_UPCXX
template<typename FAB>
FabArrayCopyDescriptor<FAB>* FabArrayCopyDescriptor<FAB>::CurrentCollectDataObj = NULL;

template<typename FAB>
std::vector<void*> FabArrayCopyDescriptor<FAB>::CollectData_send_data_buffers;

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::CollectData_AM_handler(uint32_t src_rank,
                                                    void *payload,
                                                    size_t nbytes,
                                                    void *target_buffer,
                                                    void *target_event,
                                                    int numboxes,
                                                    int isAsync)
{
  const int Nints = 4 + 3*BL_SPACEDIM;  // # of ints in a meta-data
  typedef typename FAB::value_type value_type;
  int *p = (int *)payload;

  BL_ASSERT(CurrentCollectDataObj != NULL);

  Array<int> faid(numboxes);
  Array<int> fidx(numboxes);
  Array<int> scomp(numboxes);
  Array<int> ncomp(numboxes);
  Array<int> npts(numboxes);
  Array<Box> bxs;
  int N = 0;
  const int * md = p;
  for (int i = 0; i < numboxes; ++i, md += Nints)
  {
    faid[i] = md[0];
    fidx[i] = md[1];
    scomp[i] = md[2];
    ncomp[i] = md[3];
    bxs.push_back(Box(IntVect(&md[4]),
                      IntVect(&md[4+BL_SPACEDIM]),
                      IntVect(&md[4+BL_SPACEDIM*2])));
    npts[i] = bxs.back().numPts()*ncomp[i];
    N += npts[i];
  }

  BL_ASSERT(N < std::numeric_limits<int>::max());

  value_type* data = static_cast<value_type*>(BLPgas::alloc(N*sizeof(value_type)));
  BL_ASSERT(data != NULL);
  value_type* dptr = data;
  CollectData_send_data_buffers.push_back(data);

  for (int i = 0; i < numboxes; ++i)
  {
    (*(CurrentCollectDataObj->fabArrays)[faid[i]])[fidx[i]].copyToMem(bxs[i],scomp[i],ncomp[i],dptr);
    dptr += npts[i];
  }

  if (isAsync) {
    BLPgas::free(payload);
    upcxx::async_copy_and_signal(upcxx::global_ptr<char>((char *)data),
                                 upcxx::global_ptr<char>((char *)target_buffer, src_rank),
                                 N*sizeof(value_type),
                                 (upcxx::event *)target_event,
                                 (upcxx::event *)NULL,
                                 (upcxx::event *)NULL);
  } else {
    upcxx::async(upcxx::myrank())(upcxx::async_copy_and_signal<char>,
                                  upcxx::global_ptr<char>((char *)data),
                                  upcxx::global_ptr<char>((char *)target_buffer, src_rank),
                                  N*sizeof(value_type),
                                  (upcxx::event *)target_event,
                                  (upcxx::event *)NULL,
                                  (upcxx::event *)NULL);
  }
}


template <class FAB>
void
FabArrayCopyDescriptor<FAB>::CollectData ()
{
  dataAvailable = true;

  if (ParallelDescriptor::NProcs() == 1) return;

  // The UPC++ version assumes that all processes in MPI_COMM_WORLD participate in CollectData
  BL_ASSERT(ParallelDescriptor::NProcs() == upcxx::ranks());

  typedef typename FAB::value_type value_type;
  BL_PROFILE("FabArrayCopyDescriptor::CollectData()");

  CurrentCollectDataObj = this; // store the this pointer into a static object for use in active messages
  upcxx::barrier();

  const int MyProc = ParallelDescriptor::MyProc();
  const int NProcs = ParallelDescriptor::NProcs();

  int Total_Rcvs_Size = 0;
  //
  // We use this to make finding matching FabComTags more efficient.
  //
  std::map< int,FabComTagIterContainer > RcvTags;

  std::map<int,int> Rcvs, Npts;

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

  const int N_rcvs = Rcvs.size();

  const int Nints = 4 + 3*BL_SPACEDIM;  // # of ints in a meta-data

  // for meta-data
  Array<int> md_offset, md_icnts, md_bcnts;
  int* md_recv_data;

  // for data
  Array<int> data_sender, data_offset;
  value_type* recv_data = NULL;
  upcxx::event *copy_events = NULL;

  if (N_rcvs > 0)
  {
    BLPgas::CollectData_recv_event.incref((uint32_t)N_rcvs);
    copy_events = new upcxx::event[N_rcvs];
    int event_index = 0;

    // Send meta-data
    int Idx = 0;
    recv_data = static_cast<value_type*>(BLPgas::alloc(Total_Rcvs_Size*sizeof(value_type)));
    for (std::map<int,int>::const_iterator it = Rcvs.begin(), End = Rcvs.end(); it != End; ++it)
    {
      int rank = it->first;
      int Nmds = it->second;
      int cnt = Nmds * Nints;

      int* p = static_cast<int*>(BLPgas::alloc(cnt*sizeof(int)));
      const FabComTagIterContainer& tags = RcvTags[rank];

      // initialized the data
      int * md = p;
      for (int i = 0; i < Nmds; ++i, md += Nints)
      {
        md[0] = tags[i]->fabArrayId;
        md[1] = tags[i]->fabIndex;
        md[2] = tags[i]->srcComp;
        md[3] = tags[i]->nComp;
        const int* lo = tags[i]->box.loVect();
        const int* hi = tags[i]->box.hiVect();
        const int* tp = tags[i]->box.type().getVect();
        D_EXPR(md[4] = lo[0],
               md[5] = lo[1],
               md[6] = lo[2]);
        D_EXPR(md[4+  BL_SPACEDIM] = hi[0],
               md[5+  BL_SPACEDIM] = hi[1],
               md[6+  BL_SPACEDIM] = hi[2]);
        D_EXPR(md[4+2*BL_SPACEDIM] = tp[0],
               md[5+2*BL_SPACEDIM] = tp[1],
               md[6+2*BL_SPACEDIM] = tp[2]);
      }

      data_sender.push_back(rank);
      data_offset.push_back(Idx);

      int Cnt = Npts[rank];
      BL_ASSERT(Cnt > 0);
      BL_ASSERT(Cnt < std::numeric_limits<int>::max());

      if (cnt*sizeof(int) < upcxx::max_am_payload_size()) {
        upcxx::am_send(rank, CollectData_AM_handler, p, (size_t)cnt*sizeof(int),
                       &recv_data[Idx], &BLPgas::CollectData_recv_event, Nmds, 0);
        BLPgas::free(p);
      } else
      {
        // allocate a remote buffer for the meta data because it cannot fit in an AM
        upcxx::global_ptr<int> remote_meta_data = upcxx::allocate<int>(rank, (size_t)cnt);
        BL_ASSERT(remote_meta_data != nullptr);
        upcxx::async_copy(upcxx::global_ptr<int>(p), remote_meta_data, cnt, &copy_events[event_index]);
        // execute the following two tasks after the above async_copy is done
        upcxx::async_after(rank, &copy_events[event_index])(CollectData_AM_handler,
                                                            upcxx::myrank(),
                                                            remote_meta_data.raw_ptr(),
                                                            (size_t)cnt*sizeof(int),
                                                            &recv_data[Idx],
                                                            &BLPgas::CollectData_recv_event,
                                                            Nmds, 1);
        upcxx::async_after(upcxx::myrank(), &copy_events[event_index])(BLPgas::free, p);
        event_index++;
      }
      Idx += Cnt;
    }
  }

  // Wait and unpack data
  if (N_rcvs > 0)
  {
    BLPgas::CollectData_recv_event.wait();
    std::pair<FCDMapIter,FCDMapIter> match;
    std::map< int,FabComTagIterContainer >::const_iterator found;
    for (int k = 0; k < N_rcvs; k++)
    {
      const int         Who  = data_sender[k];
      const value_type* dptr = &recv_data[data_offset[k]];
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
    BLPgas::free(recv_data);
  }
  upcxx::async_wait();
  if (copy_events) delete[] copy_events;

  upcxx::barrier();
  ParallelDescriptor::Mode.set_mpi_mode(true); // since we just had a barrier

  // It's safe to free the send buffers now
  for (auto it = CollectData_send_data_buffers.begin(); it != CollectData_send_data_buffers.end(); ++it) {
    if (*it != NULL) BLPgas::free(*it);
  }
  CollectData_send_data_buffers.clear();
  CurrentCollectDataObj = NULL;
}

#else // else of #if BL_USE_UPCXX


template <class FAB>
void
FabArrayCopyDescriptor<FAB>::CollectData ()
{
    dataAvailable = true;

    if (ParallelDescriptor::NProcs() == 1) return;

#if BL_USE_MPI
    typedef typename FAB::value_type value_type;

    BL_PROFILE("FabArrayCopyDescriptor::CollectData()");

    const int MyProc = ParallelDescriptor::MyProc();

    int Total_Rcvs_Size = 0;
    //
    // We use this to make finding matching FabComTags more efficient.
    //
    std::map< int,FabComTagIterContainer > RcvTags;

    std::map<int,int> Snds, Rcvs, Npts;
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

        for (std::map<int,int>::const_iterator it = Rcvs.begin(), End = Rcvs.end(); it != End; ++it)
	{
            RcvsArray[it->first] = it->second;
	}

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

    // There are two rounds of send and recv.
    // First, the data receivers need to send the data senders meta-data (e.g., boxes).
    // Then, the senders know what data to send and perform send.
    const int SeqNum_md   = ParallelDescriptor::SeqNum();
    const int SeqNum_data = ParallelDescriptor::SeqNum();

    const int N_snds = Snds.size();
    const int N_rcvs = Rcvs.size();

    if ( N_snds == 0 && N_rcvs == 0 ) return;

    const int Nints = 4 + 3*BL_SPACEDIM;  // # of ints in a meta-data

    // for meta-data
    Array<int> md_sender, md_offset, md_icnts, md_bcnts;
    int* md_recv_data;
    Array<int*> md_send_data;
    Array<MPI_Request> md_recv_reqs, md_send_reqs;

    // for data
    Array<int> data_sender, data_offset;
    value_type* recv_data;
    Array<value_type*> send_data;
    Array<MPI_Request> data_recv_reqs, data_send_reqs;

    if (N_snds > 0)
    {
	// Recv meta-data

	int N = 0;
        for (std::map<int,int>::const_iterator it = Snds.begin(), End = Snds.end(); it != End; ++it)
        {
	    md_sender.push_back(it->first);
	    md_bcnts.push_back(it->second);
	    int cnt = it->second * Nints;
	    md_icnts.push_back(cnt);
            md_offset.push_back(N);
            N += cnt; 
        }

	md_recv_data = static_cast<int*>(BoxLib::The_Arena()->alloc(N*sizeof(int)));

	for (int i = 0; i < N_snds; ++i)
	{
	    md_recv_reqs.push_back(ParallelDescriptor::Arecv(&md_recv_data[md_offset[i]],
							     md_icnts[i], md_sender[i],
							     SeqNum_md).req());
	}	
    }

    if (N_rcvs > 0)
    {
	// Send meta-data
        for (std::map<int,int>::const_iterator it = Rcvs.begin(), End = Rcvs.end(); it != End; ++it)
        {
	    int rank = it->first;
	    int Nmds = it->second;
	    int cnt = Nmds * Nints;

	    int* p = static_cast<int*>(BoxLib::The_Arena()->alloc(cnt*sizeof(int)));
	    md_send_data.push_back(p);

	    const FabComTagIterContainer& tags = RcvTags[rank];

	    // initialized the data
	    int * md = p;
	    for (int i = 0; i < Nmds; ++i, md += Nints)
	    {
		md[0] = tags[i]->fabArrayId;
		md[1] = tags[i]->fabIndex;
		md[2] = tags[i]->srcComp;
		md[3] = tags[i]->nComp;
		const int* lo = tags[i]->box.loVect();
		const int* hi = tags[i]->box.hiVect();
		const int* tp = tags[i]->box.type().getVect();
		D_EXPR(md[4] = lo[0],
		       md[5] = lo[1],
		       md[6] = lo[2]);
		D_EXPR(md[4+  BL_SPACEDIM] = hi[0],
		       md[5+  BL_SPACEDIM] = hi[1],
		       md[6+  BL_SPACEDIM] = hi[2]);
		D_EXPR(md[4+2*BL_SPACEDIM] = tp[0],
		       md[5+2*BL_SPACEDIM] = tp[1],
		       md[6+2*BL_SPACEDIM] = tp[2]);
	    }

	    md_send_reqs.push_back(ParallelDescriptor::Asend(p,cnt,rank,SeqNum_md).req());
	}
    }

    if (N_rcvs > 0)
    {
	recv_data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(Total_Rcvs_Size*sizeof(value_type)));

	// Post receives for data
	int Idx = 0;
	for (std::map<int,int>::const_iterator it = Npts.begin(); it != Npts.end(); ++it)
	{
	    int Who = it->first;
	    int Cnt = it->second;
	    BL_ASSERT(Cnt > 0);
	    BL_ASSERT(Cnt < std::numeric_limits<int>::max());
	    data_sender.push_back(Who);
	    data_recv_reqs.push_back(ParallelDescriptor::Arecv(&recv_data[Idx],
							       Cnt,Who,SeqNum_data).req());
	    data_offset.push_back(Idx);
	    Idx += Cnt;
	}
    }

    // Wait on meta-data and do send
    if (N_snds > 0) 
    {
	int send_counter = 0;
	while (send_counter++ < N_snds)
	{
	    MPI_Status status;
	    int index;
	    MPI_Waitany(N_snds, md_recv_reqs.dataPtr(), &index, &status);

	    int rank = status.MPI_SOURCE;
	    BL_ASSERT(status.MPI_TAG == SeqNum_md);
	    BL_ASSERT(rank == md_sender[index]);

	    const int* p = &md_recv_data[md_offset[index]];
	    int numboxes = md_bcnts[index];
	    Array<int> faid(numboxes);
	    Array<int> fidx(numboxes);
	    Array<int> scomp(numboxes);
	    Array<int> ncomp(numboxes);
	    Array<int> npts(numboxes);
	    Array<Box> bxs;
	    int N = 0;
	    const int * md = p;
	    for (int i = 0; i < numboxes; ++i, md += Nints)
	    {
		faid[i] = md[0];
		fidx[i] = md[1];
		scomp[i] = md[2];
		ncomp[i] = md[3];
		bxs.push_back(Box(IntVect(&md[4]),
				  IntVect(&md[4+BL_SPACEDIM]), 
				  IntVect(&md[4+BL_SPACEDIM*2]))); 
		npts[i] = bxs.back().numPts()*ncomp[i];
		N += npts[i];
	    }

	    BL_ASSERT(N < std::numeric_limits<int>::max());
	    
	    value_type* data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(N*sizeof(value_type)));
	    value_type* dptr = data;
	    send_data.push_back(data);

	    for (int i = 0; i < numboxes; ++i)
	    {
		(*fabArrays[faid[i]])[fidx[i]].copyToMem(bxs[i],scomp[i],ncomp[i],dptr);
		dptr += npts[i];
	    }	    

	    data_send_reqs.push_back(ParallelDescriptor::Asend(data,N,rank,SeqNum_data).req());
	}

	BoxLib::The_Arena()->free(md_recv_data);
    }

    // Wait and upack data
    if (N_rcvs > 0)
    {
	Array<MPI_Status> stats(N_rcvs);

	BL_MPI_REQUIRE (MPI_Waitall(N_rcvs, md_send_reqs.dataPtr(), stats.dataPtr()) );
	for (int i = 0; i < N_rcvs; ++i) {
            BoxLib::The_Arena()->free(md_send_data[i]);	    
	}

	BL_MPI_REQUIRE( MPI_Waitall(N_rcvs, data_recv_reqs.dataPtr(), stats.dataPtr()) );

	std::pair<FCDMapIter,FCDMapIter> match;
	std::map< int,FabComTagIterContainer >::const_iterator found;

	for (int k = 0; k < N_rcvs; k++)
	{
	    const int         Who  = data_sender[k];
	    const value_type* dptr = &recv_data[data_offset[k]];
	    
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

	BoxLib::The_Arena()->free(recv_data);
    }

    // Finished send
    if (N_snds > 0)
    {
	Array<MPI_Status> stats(N_snds);
        BL_COMM_PROFILE(BLProfiler::Waitall, sizeof(value_type), BLProfiler::BeforeCall(), N_snds);
	BL_MPI_REQUIRE( MPI_Waitall(N_snds, data_send_reqs.dataPtr(), stats.dataPtr()) );
        BL_COMM_PROFILE(BLProfiler::Waitall, sizeof(value_type), BLProfiler::AfterCall(), N_snds);

	for (int i = 0; i < N_snds; ++i) {
            BoxLib::The_Arena()->free(send_data[i]);
	}
    }

#endif /*BL_USE_MPI*/
}
#endif // end of !BL_USE_UPCXX

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
#pragma omp parallel for if (TheFB.m_threadsafe_loc)
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

#if !defined(BL_USE_MPI3)
    BL_ASSERT(!ParallelDescriptor::MPIOneSided());
#endif

    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    int SeqNum;
    {
	ParallelDescriptor::Color mycolor = this->color();
	if (mycolor == ParallelDescriptor::DefaultColor()) {
	    SeqNum = ParallelDescriptor::SeqNum();
	} else if (mycolor == ParallelDescriptor::SubCommColor()) {
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
    // Post rcvs. Allocate one chunk of space to hold'm all.
    //
#ifdef BL_USE_MPI3
    MPI_Group tgroup, rgroup, sgroup;
    if (ParallelDescriptor::MPIOneSided())
	MPI_Comm_group(ParallelDescriptor::Communicator(), &tgroup);
#endif

    if (N_rcvs > 0) {
#ifdef BL_USE_UPCXX
	FabArrayBase::PostRcvs_PGAS(*TheFB.m_RcvVols,fb_the_recv_data,
				    fb_recv_data,fb_recv_from,ncomp,SeqNum,&BLPgas::fb_recv_event);
#else
	if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
	    FabArrayBase::PostRcvs_MPI_Onesided(*TheFB.m_RcvVols, fb_the_recv_data,
						fb_recv_data, fb_recv_from, fb_recv_reqs, 
						fb_recv_disp, ncomp, SeqNum, ParallelDescriptor::fb_win);
	    MPI_Group_incl(tgroup, fb_recv_from.size(), fb_recv_from.dataPtr(), &rgroup);
	    MPI_Win_post(rgroup, 0, ParallelDescriptor::fb_win);
#endif
	} else {
	    FabArrayBase::PostRcvs(*TheFB.m_RcvVols,fb_the_recv_data,
				   fb_recv_data,fb_recv_from,fb_recv_reqs,ncomp,SeqNum);
	}
#endif
    }

    //
    // Post send's
    //
    if (N_snds > 0)
    {
        Array<value_type*> &               send_data = fb_send_data;
	Array<int>                         send_N;
	Array<int>                         send_rank;
	Array<const CopyComTagsContainer*> send_cctc;

	send_data.reserve(N_snds);
	send_N   .reserve(N_snds);
	send_rank.reserve(N_snds);
	send_cctc.reserve(N_snds);

	for (MapOfCopyComTagContainers::const_iterator m_it = TheFB.m_SndTags->begin(),
		 m_End = TheFB.m_SndTags->end();
	     m_it != m_End;
	     ++m_it)
	{
	    std::map<int,int>::const_iterator vol_it = TheFB.m_SndVols->find(m_it->first);

	    BL_ASSERT(vol_it != TheFB.m_SndVols->end());
	    
	    const int N = vol_it->second*ncomp;
	    
	    BL_ASSERT(N < std::numeric_limits<int>::max());
	    
	    value_type* data = static_cast<value_type*>
#ifdef BL_USE_UPCXX
		(BLPgas::alloc(N*sizeof(value_type)));
#else
	        (BoxLib::The_Arena()->alloc(N*sizeof(value_type)));
#endif

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
		const Box& bx = it->sbox;
		get(it->srcIndex).copyToMem(bx,scomp,ncomp,dptr);
		const int Cnt = bx.numPts()*ncomp;
		dptr += Cnt;
	    }
	}

#ifdef BL_USE_UPCXX

	BLPgas::fb_send_counter = 0;

	for (int i=0; i<N_snds; ++i) {
	    BLPgas::Send(upcxx::global_ptr<void>((void *)send_data[i], upcxx::myrank()),
			 send_rank[i],
			 send_N[i]*sizeof(value_type),
			 SeqNum,
			 &BLPgas::fb_send_event,
			 &BLPgas::fb_send_counter);
	}

	// Need to make sure at least half of the sends have been started
	while (BLPgas::fb_send_counter < N_snds)
	    upcxx::advance();

#else  // MPI

	if (ParallelDescriptor::MPIOneSided())
	{
#if defined(BL_USE_MPI3)
	    Array<MPI_Request> send_reqs;
	    Array<MPI_Aint>    send_disp;
	    send_reqs.reserve(N_snds);
	    send_disp.resize (N_snds);
	    
	    for (int i=0; i<N_snds; ++i) {
		send_reqs.push_back(ParallelDescriptor::Arecv
				    (&send_disp[i],1,send_rank[i],SeqNum).req());
	    }
		
	    MPI_Group_incl(tgroup, send_rank.size(), send_rank.dataPtr(), &sgroup);
	    MPI_Win_start(sgroup,0,ParallelDescriptor::fb_win);
	    
	    int send_counter = 0;	
	    while (send_counter < N_snds) {
		MPI_Status status;
		int index;
		
		MPI_Waitany(N_snds, send_reqs.dataPtr(), &index, &status);
		
		BL_ASSERT(status.MPI_TAG == SeqNum);
		BL_ASSERT(status.MPI_SOURCE == send_rank[index]);
		
		MPI_Put(send_data[index], send_N[index]*sizeof(value_type), MPI_CHAR, send_rank[index],
			send_disp[index], send_N[index]*sizeof(value_type), MPI_CHAR, 
			ParallelDescriptor::fb_win);
		
		send_counter++;
	    }
#endif // BL_USE_MPI3
	} 
	else 
	{
	    fb_send_reqs.reserve(N_snds);

	    for (int i=0; i<N_snds; ++i) {
		fb_send_reqs.push_back(ParallelDescriptor::Asend
				       (send_data[i],send_N[i],send_rank[i],SeqNum).req());
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
#pragma omp parallel
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
#pragma omp parallel for if (TheFB.m_threadsafe_loc)
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

#ifdef BL_USE_UPCXX
    if (N_rcvs > 0) BLPgas::fb_recv_event.wait();
#else 
    if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
	if (N_snds > 0) MPI_Win_complete(ParallelDescriptor::fb_win);
	if (N_rcvs > 0) MPI_Win_wait    (ParallelDescriptor::fb_win);
#endif
    } else {
	if (N_rcvs > 0) {
	    Array<MPI_Status> stats(N_rcvs);
	    BL_MPI_REQUIRE( MPI_Waitall(N_rcvs, fb_recv_reqs.dataPtr(), stats.dataPtr()) );
	}
    }
#endif

    if (N_rcvs > 0)
    {
	Array<const CopyComTagsContainer*> recv_cctc;
	recv_cctc.reserve(N_rcvs);

	for (int k = 0; k < N_rcvs; k++) 
	{
	    MapOfCopyComTagContainers::const_iterator m_it = TheFB.m_RcvTags->find(fb_recv_from[k]);
	    BL_ASSERT(m_it != TheFB.m_RcvTags->end());
	    recv_cctc.push_back(&(m_it->second));
	}	

#ifdef _OPENMP
#pragma omp parallel for if (TheFB.m_threadsafe_rcv)
#endif
	for (int k = 0; k < N_rcvs; k++) 
	{
	    value_type* dptr = fb_recv_data[k];
	    BL_ASSERT(dptr != 0);

	    const CopyComTagsContainer& cctc = *recv_cctc[k];
	    
	    for (CopyComTagsContainer::const_iterator it = cctc.begin();
		 it != cctc.end(); ++it)
	    {
		const Box& bx  = it->dbox;
		const int  Cnt = bx.numPts()*fb_ncomp;
		get(it->dstIndex).copyFromMem(bx,fb_scomp,fb_ncomp,dptr);
		dptr += Cnt;
	    }	    
	}

#ifdef BL_USE_UPCXX
	BLPgas::free(fb_the_recv_data);
#else
	if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
	    MPI_Win_detach(ParallelDescriptor::fb_win, fb_the_recv_data);
	    BoxLib::The_Arena()->free(fb_the_recv_data);
	    fb_recv_disp.clear();
#endif
	} else {
	    BoxLib::The_Arena()->free(fb_the_recv_data);
	}
#endif

	fb_recv_from.clear();
	fb_recv_data.clear();
	fb_recv_reqs.clear();
    }

    if (N_snds > 0) {
#ifdef BL_USE_UPCXX
        FabArrayBase::WaitForAsyncSends_PGAS(N_snds,fb_send_data,
					     &BLPgas::fb_send_event,
					     &BLPgas::fb_send_counter);
#else
	if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
	    for (int i = 0; i < N_snds; ++i)
		BoxLib::The_Arena()->free(fb_send_data[i]);
#endif
	} else {
	    Array<MPI_Status> stats;
	    FabArrayBase::WaitForAsyncSends(N_snds,fb_send_reqs,fb_send_data,stats);
	}
#endif
	fb_send_data.clear();
	fb_send_reqs.clear();
    }

#ifdef BL_USE_TEAM
    ParallelDescriptor::MyTeam().MemoryBarrier();
#endif

#endif // MPI
}

#ifdef BL_USE_UPCXX
template<typename T>
void
FabArrayBase::PostRcvs_PGAS (const std::map<int,int>&               m_RcvVols,
                             T*&                                    the_recv_data,
                             Array<T*>&                             recv_data,
                             Array<int>&                            recv_from,
                             int                                    ncomp,
                             int                                    SeqNum,
                             upcxx::event*                          recv_event)
{
  int TotalRcvsVolume = 0;

  for (auto it = m_RcvVols.begin(); it != m_RcvVols.end(); ++it)
  {
      TotalRcvsVolume += it->second;
  }

  TotalRcvsVolume *= ncomp;

  BL_ASSERT((TotalRcvsVolume*sizeof(T)) < std::numeric_limits<int>::max());

  the_recv_data = (T*) BLPgas::alloc(TotalRcvsVolume*sizeof(T));

  int Offset = 0;

  for (auto it = m_RcvVols.begin(); it != m_RcvVols.end(); ++it)
  {
      const int N = it->second*ncomp;

      BL_ASSERT(N < std::numeric_limits<int>::max());

      recv_data.push_back(&the_recv_data[Offset]);
      recv_from.push_back(it->first);
      /// Send an AM to the sender it->first with the recv pointer being &the_recv_data[Offset]
      upcxx::global_ptr<void> dst_ptr =
        upcxx::global_ptr<void>(&the_recv_data[Offset], upcxx::myrank());
      // Increment the event reference before launching the remote task
      recv_event->incref();
      // Launch a remote task on the sender to send me data
      BLPgas::Request(it->first,
                      dst_ptr,
                      N*sizeof(T),
                      SeqNum,
                      recv_event);
      upcxx::advance(); // poll the UPC++ progress engine and the network
      Offset += N;
    }
}

template<typename T>
void
FabArrayBase::WaitForAsyncSends_PGAS (int                 N_snds,
                                      Array<T*>&          send_data,
                                      upcxx::event*       send_event,
                                      volatile int*       send_counter)
{
  BL_ASSERT(N_snds > 0);
  BL_ASSERT(send_data.size() == N_snds);
  // Need to make sure all sends have been started
  while ((*send_counter) < N_snds)
      upcxx::advance();
  send_event->wait(); // wait for the sends
  for (int i = 0; i < N_snds; i++)
      BLPgas::free(send_data[i]);
}
#endif /* BL_USE_UPCXX */


template <class FAB>
void
FabArray<FAB>::BuildMask (const Box& phys_domain, const Periodicity& period,
			  value_type covered, value_type uncovered, 
			  value_type physbnd, value_type interior)
{
    int ncomp = this->nComp();
    int ngrow = this->nGrow();

    const FabArrayBase::FB& TheFB = this->getFB(period);

    const CopyComTagsContainer&      LocTags = *(TheFB.m_LocTags);
    const MapOfCopyComTagContainers& RcvTags = *(TheFB.m_RcvTags);

    Box domain = BoxLib::convert(phys_domain, boxArray().ixType());
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	if (period.isPeriodic(i)) {
	    domain.grow(i, ngrow);
	}
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
	FAB& fab = (*this)[mfi];

	Box gbx = mfi.growntilebox();
	fab.setVal(physbnd, gbx, 0, ncomp);

	gbx &= domain;
	fab.setVal(uncovered, gbx, 0, ncomp);

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

#endif /*BL_FABARRAY_H*/
