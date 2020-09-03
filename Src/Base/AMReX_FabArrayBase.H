#ifndef BL_FABARRAYBASE_H_
#define BL_FABARRAYBASE_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <string>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Periodicity.H>
#include <AMReX_Print.H>
#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>

namespace amrex {

class MFIter;
class MFGhostIter;
class Geometry;
class FArrayBox;
template <typename FAB> class FabFactory;
template <typename FAB> class FabArray;

namespace EB2 { class IndexSpace; }

class FabArrayBase
{
    friend class MFIter;
    friend class MFGhostIter;

    template <class FAB> friend void FillBoundary (Vector<FabArray<FAB>*> const& mf, const Periodicity& period);

public:

    FabArrayBase ();

    FabArrayBase (const BoxArray&            bxs,
                  const DistributionMapping& dm,
                  int                        nvar,
                  int                        ngrow);

    FabArrayBase (const BoxArray&            bxs,
                  const DistributionMapping& dm,
                  int                        nvar,
                  const IntVect&             ngrow);

    virtual ~FabArrayBase();

    FabArrayBase (FabArrayBase&& rhs) noexcept = default;
    FabArrayBase (const FabArrayBase& rhs) = default;
    FabArrayBase& operator= (const FabArrayBase& rhs) = default;
    FabArrayBase& operator= (FabArrayBase&& rhs) = default;

    void define (const BoxArray&            bxs,
                 const DistributionMapping& dm,
                 int                        nvar,
                 int                        ngrow);

    void define (const BoxArray&            bxs,
                 const DistributionMapping& dm,
                 int                        nvar,
                 const IntVect&             ngrow);

    //! Return the grow factor that defines the region of definition.
    int nGrow (int direction = 0) const noexcept { return n_grow[direction]; }

    IntVect nGrowVect () const noexcept { return n_grow; }

    //! Return number of variables (aka components) associated with each point.
    int nComp () const noexcept { return n_comp; }

    //! Return index type.
    IndexType ixType () const noexcept { return boxarray.ixType(); }

    //Return whether this FabArray is empty
    bool empty () const noexcept { return boxarray.empty(); }

    /**
    * \brief Return a constant reference to the BoxArray that defines the
    * valid region associated with this FabArray.
    */
    const BoxArray& boxArray () const noexcept { return boxarray; }

    /**
    * \brief Return the Kth Box in the BoxArray.
    * That is, the valid region of the Kth grid.
    */
    Box box (int K) const noexcept { return boxarray[K]; }

    /**
    * \brief Return the Kth FABs Box in the FabArray.
    * That is, the region the Kth fab is actually defined on.
    */
    Box fabbox (int K) const noexcept;

    //! Return the number of FABs in the FabArray.
    int size () const noexcept { return boxarray.size(); }

    //! Return the number of local FABs in the FabArray.
    int local_size () const noexcept { return indexArray.size(); }

    //! Return constant reference to indices in the FabArray that we have access.
    const Vector<int> &IndexArray () const noexcept { return indexArray; }

    //! Return local index in the vector of FABs.
    int localindex (int K) const noexcept {
        std::vector<int>::const_iterator low
            = std::lower_bound(indexArray.begin(), indexArray.end(), K);
        if (low != indexArray.end() && *low == K) {
            return low - indexArray.begin();
        }
        else {
            return -1;
        }
    }

    //! Return constant reference to associated DistributionMapping.
    const DistributionMapping& DistributionMap () const noexcept { return distributionMap; }

    /**
    * \brief This tests on whether the FabArray is fully nodal.
    */
    bool is_nodal () const noexcept;
    /**
    * \brief This tests on whether the FabArray is nodal in direction dir.
    */
    bool is_nodal (int dir) const noexcept;
    /**
    * \brief This tests on whether the FabArray is cell-centered.
    */
    bool is_cell_centered () const noexcept;

    void setMultiGhost(bool a_multi_ghost) {m_multi_ghost = a_multi_ghost;}

    // These are provided for convenience to keep track of how many
    // ghost cells are up to date.  The number of filled ghost cells
    // is updated by FillBoundary and ParallelCopy.
    IntVect nGrowFilled () const noexcept { return n_filled; }
    void setNGrowFilled (IntVect const& ng) noexcept { n_filled = ng; }

    //
    struct CacheStats
    {
	int         size;     //!< current size: nbuild - nerase
	int         maxsize;  //!< highest water mark of size
	Long        maxuse;   //!< max # of uses of a cached item
	Long        nuse;     //!< # of uses of the whole cache
	Long        nbuild;   //!< # of build operations
	Long        nerase;   //!< # of erase operations
	Long        bytes;
	Long        bytes_hwm;
	std::string name;     //!< name of the cache
	explicit CacheStats (const std::string& name_)
	    : size(0),maxsize(0),maxuse(0),nuse(0),nbuild(0),nerase(0),
	      bytes(0L),bytes_hwm(0L),name(name_) {;}
	void recordBuild () noexcept {
	    ++size;
	    ++nbuild;
	    maxsize = std::max(maxsize, size);
	}
	void recordErase (Long n) noexcept {
	    // n: how many times the item to be deleted has been used.
	    --size;
	    ++nerase;
	    maxuse = std::max(maxuse, n);
	}
	void recordUse () noexcept { ++nuse; }
	void print () {
	    amrex::Print(Print::AllProcs) << "### " << name << " ###\n"
					  << "    tot # of builds  : " << nbuild  << "\n"
					  << "    tot # of erasures: " << nerase  << "\n"
					  << "    tot # of uses    : " << nuse    << "\n"
					  << "    max cache size   : " << maxsize << "\n"
					  << "    max # of uses    : " << maxuse  << "\n";
	}
    };
    //
    //! Used by a bunch of routines when communicating via MPI.
    struct CopyComTag
    {
	Box dbox;
        Box sbox;
        int dstIndex;
        int srcIndex;
	CopyComTag () noexcept {}
	CopyComTag (const Box& db, const Box& sb, int didx, int sidx) noexcept
	    : dbox(db), sbox(sb), dstIndex(didx), srcIndex(sidx) {}
	bool operator< (const CopyComTag& rhs) const noexcept {
	    return (srcIndex < rhs.srcIndex) || ((srcIndex == rhs.srcIndex) && (
                   (sbox.smallEnd() < rhs.sbox.smallEnd()
		               || ((sbox.smallEnd() == rhs.sbox.smallEnd()) && (
                   (dstIndex < rhs.dstIndex) || ((dstIndex == rhs.dstIndex) && (
                   (dbox.smallEnd() < rhs.dbox.smallEnd()))))))));
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
    static Long bytesOfMapOfCopyComTagContainers (const MapOfCopyComTagContainers&);

    /**
    * Key for unique combination of BoxArray and DistributionMapping
    * Note both BoxArray and DistributionMapping are reference counted.
    * Objects with the same references have the same key.
    */
    struct BDKey {
        BDKey () noexcept = default;
        BDKey (const BoxArray::RefID& baid, const DistributionMapping::RefID& dmid) noexcept
            : m_ba_id(baid), m_dm_id(dmid) {}
        bool operator< (const BDKey& rhs) const noexcept {
            return (m_ba_id < rhs.m_ba_id) ||
                ((m_ba_id == rhs.m_ba_id) && (m_dm_id < rhs.m_dm_id));
        }
        bool operator== (const BDKey& rhs) const noexcept {
            return m_ba_id == rhs.m_ba_id && m_dm_id == rhs.m_dm_id;
        }
        bool operator!= (const BDKey& rhs) const noexcept {
            return m_ba_id != rhs.m_ba_id || m_dm_id != rhs.m_dm_id;
        }
        friend std::ostream& operator<< (std::ostream& os, const BDKey& id);
    private:
        BoxArray::RefID            m_ba_id;
        DistributionMapping::RefID m_dm_id;
    };

    BDKey getBDKey () const noexcept {
	return {boxarray.getRefID(), distributionMap.getRefID()};
    }

    void updateBDKey ();

    //
    //! Tiling
    struct TileArray
    {
        Long nuse;
        Vector<int> numLocalTiles;
        Vector<int> indexMap;
        Vector<int> localIndexMap;
        Vector<int> localTileIndexMap;
        Vector<Box> tileArray;
        TileArray () noexcept : nuse(-1) {;}
        Long bytes () const;
    };

    //
    //! Used for collecting information used in communicating FABs.
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

        FabComTag () noexcept
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

    //! Default tilesize in MFIter
    static IntVect mfiter_tile_size;

    //! Default tilesize in MFGhostIter
    static IntVect mfghostiter_tile_size;

    //! The maximum number of components to copy() at a time.
    static int MaxComp;

    //! Initialize from ParmParse with "fabarray" prefix.
    static void Initialize ();
    static void Finalize ();
    /**
    * To maximize thread efficiency we now can decompose things like
    * intersections among boxes into smaller tiles. This sets
    * their maximum size.
    */
    static IntVect comm_tile_size;  //!< communication tile size

    struct FPinfo
    {
        FPinfo (const FabArrayBase& srcfa,
                const FabArrayBase& dstfa,
                const Box&          dstdomain,
                const IntVect&      dstng,
                const BoxConverter& coarsener,
                const Box&          fdomain,
                const Box&          cdomain,
                const EB2::IndexSpace* index_space);

        ~FPinfo ();

        Long bytes () const;

        BoxArray            ba_crse_patch;
        BoxArray            ba_fine_patch;
        DistributionMapping dm_patch;
        std::unique_ptr<FabFactory<FArrayBox> > fact_crse_patch;
        std::unique_ptr<FabFactory<FArrayBox> > fact_fine_patch;
        //
        BDKey               m_srcbdk;
        BDKey               m_dstbdk;
        Box                 m_dstdomain;
        IntVect             m_dstng;
        std::unique_ptr<BoxConverter> m_coarsener;
        //
        Long                m_nuse;
    };

    typedef std::multimap<BDKey,FabArrayBase::FPinfo*> FPinfoCache;
    typedef FPinfoCache::iterator FPinfoCacheIter;

    static FPinfoCache m_TheFillPatchCache;

    static CacheStats m_FPinfo_stats;

    static const FPinfo& TheFPinfo (const FabArrayBase& srcfa,
                                    const FabArrayBase& dstfa,
                                    const IntVect&      dstng,
                                    const BoxConverter& coarsener,
                                    const Geometry&     fgeom,
                                    const Geometry&     cgeom,
                                    const EB2::IndexSpace*);

    void flushFPinfo (bool no_assertion=false);

    //
    //! coarse/fine boundary
    struct CFinfo
    {
        CFinfo (const FabArrayBase& finefa,
                const Geometry&     finegm,
                const IntVect&      ng,
                bool                include_periodic,
                bool                include_physbndry);

        Long bytes () const;

        static Box Domain (const Geometry& geom, const IntVect& ng,
                           bool include_periodic, bool include_physbndry);

        BoxArray            ba_cfb;
        DistributionMapping dm_cfb;
        Vector<int>          fine_grid_idx; //!< local array
        //
        BDKey               m_fine_bdk;
        Box                 m_fine_domain;
        IntVect             m_ng;
        bool                m_include_periodic;
        bool                m_include_physbndry;
        //
        Long                m_nuse;
    };

    using CFinfoCache = std::multimap<BDKey,FabArrayBase::CFinfo*>;
    using CFinfoCacheIter = CFinfoCache::iterator;

    static CFinfoCache m_TheCrseFineCache;

    static CacheStats m_CFinfo_stats;

    static const CFinfo& TheCFinfo (const FabArrayBase& finefa,
                                    const Geometry&     finegm,
                                    const IntVect&      ng,
                                    bool                include_periodic,
                                    bool                include_physbndry);

    void flushCFinfo (bool no_assertion=false);

    //
    //! parallel copy or add
    enum CpOp { COPY = 0, ADD = 1 };

    const TileArray* getTileArray (const IntVect& tilesize) const;

    //! Block until all send requests complete
    static void WaitForAsyncSends (int                 N_snds,
                                   Vector<MPI_Request>& send_reqs,
                                   Vector<char*>&       send_data,
                                   Vector<MPI_Status>&  stats);

    // Memory Usage Tags
    struct meminfo {
        Long nbytes = 0L;
        Long nbytes_hwm = 0L;
    };
    static std::map<std::string, meminfo> m_mem_usage;

    static void updateMemUsage (std::string const& tag, Long nbytes, Arena const* ar);
    static void printMemUsage ();
    static Long queryMemUsage (const std::string& tag = std::string("All"));
    static Long queryMemUsageHWM (const std::string& tag = std::string("All"));

    static void pushRegionTag (const char* t);
    static void pushRegionTag (std::string t);
    static void popRegionTag ();

    static std::vector<std::string> m_region_tag;
    struct RegionTag {
        RegionTag (const char* t) { pushRegionTag(t); }
        RegionTag (std::string t) { pushRegionTag(std::move(t)); }
        ~RegionTag () { popRegionTag(); }
    };

//#ifndef AMREX_USE_GPU
//protected:
//#endif

    void clear ();

    /**
    * \brief Return owenership of fabs. The concept of ownership only applies when UPC++
    * team is used. In that case, each fab is shared by team workers, with one
    * taking the ownership.
    */
    const std::vector<bool>& OwnerShip () const noexcept { return ownership; }
    bool isOwner (int li) const noexcept { return ownership[li]; }

    //
    // The data ...
    //
    mutable BoxArray    boxarray;
    DistributionMapping distributionMap;
    Vector<int>         indexArray;
    std::vector<bool>   ownership;
    IntVect             n_grow;
    int                 n_comp;
    mutable BDKey       m_bdkey;
    IntVect             n_filled;  // Note that IntVect is zero by default.
    bool                m_multi_ghost = false;

    //
    // Tiling
    //
    // We use tile size as the key for the inner map.

    using TAMap   = std::map<std::pair<IntVect,IntVect>, TileArray>;
    using TACache = std::map<BDKey, TAMap>;
    //
    static TACache     m_TheTileArrayCache;
    static CacheStats  m_TAC_stats;
    //
    void buildTileArray (const IntVect& tilesize, TileArray& ta) const;
    //
    void flushTileArray (const IntVect& tilesize = IntVect::TheZeroVector(),
			 bool no_assertion=false) const;
    static void flushTileArrayCache (); //!< This flushes the entire cache.

    struct CommMetaData
    {
        // The cache of local and send/recv per FillBoundary() or ParallelCopy().
	bool m_threadsafe_loc = false;
	bool m_threadsafe_rcv = false;
        std::unique_ptr<CopyComTagsContainer>      m_LocTags;
        std::unique_ptr<MapOfCopyComTagContainers> m_SndTags;
        std::unique_ptr<MapOfCopyComTagContainers> m_RcvTags;
    };

    //
    //! FillBoundary
    struct FB
        : CommMetaData
    {
        FB (const FabArrayBase& fa, const IntVect& nghost,
            bool cross, const Periodicity& period,
	    bool enforce_periodicity_only, bool multi_ghost = false);
        ~FB ();

        IndexType    m_typ;
        IntVect      m_crse_ratio; //!< BoxArray in FabArrayBase may have crse_ratio.
        IntVect      m_ngrow;
        bool         m_cross;
        bool         m_epo;
        Periodicity  m_period;
        //
        Long         m_nuse;
        bool         m_multi_ghost = false;
        //
#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10) )
        CudaGraph<CopyMemory> m_localCopy;
        CudaGraph<CopyMemory> m_copyToBuffer;
        CudaGraph<CopyMemory> m_copyFromBuffer;
#endif
        //
        Long bytes () const;
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
    const FB& getFB (const IntVect& nghost, const Periodicity& period,
                     bool cross=false, bool enforce_periodicity_only = false) const;
    //
    void flushFB (bool no_assertion=false) const;       //!< This flushes its own FB.
    static void flushFBCache (); //!< This flushes the entire cache.

    //
    //! parallel copy or add
    struct CPC
        : CommMetaData
    {
        CPC (const FabArrayBase& dstfa, const IntVect& dstng,
             const FabArrayBase& srcfa, const IntVect& srcng,
             const Periodicity& period);
        CPC (const BoxArray& dstba, const DistributionMapping& dstdm,
             const Vector<int>& dstidx, const IntVect& dstng,
             const BoxArray& srcba, const DistributionMapping& srcdm,
             const Vector<int>& srcidx, const IntVect& srcng,
             const Periodicity& period, int myproc);
        CPC (const BoxArray& ba, const IntVect& ng,
             const DistributionMapping& dstdm, const DistributionMapping& srcdm);
        ~CPC ();

        Long bytes () const;

        BDKey       m_srcbdk;
        BDKey       m_dstbdk;
        IntVect     m_srcng;
        IntVect     m_dstng;
        Periodicity m_period;
        BoxArray    m_srcba;
        BoxArray    m_dstba;
        //
        Long        m_nuse;

    private:
        void define (const BoxArray& ba_dst, const DistributionMapping& dm_dst,
                     const Vector<int>& imap_dst,
                     const BoxArray& ba_src, const DistributionMapping& dm_src,
                     const Vector<int>& imap_src,
                     int MyProc = ParallelDescriptor::MyProc());
    };

    //
    typedef std::multimap<BDKey,FabArrayBase::CPC*> CPCache;
    typedef CPCache::iterator CPCacheIter;
    //
    static CPCache    m_TheCPCache;
    static CacheStats m_CPC_stats;
    //
    const CPC& getCPC (const IntVect& dstng, const FabArrayBase& src, const IntVect& srcng,
                       const Periodicity& period) const;
    //
    void flushCPC (bool no_assertion=false) const;      //!< This flushes its own CPC.
    static void flushCPCache (); //!< This flusheds the entire cache.

    //
    //! Keep track of how many FabArrays are built with the same BDKey.
    static std::map<BDKey, int> m_BD_count;
    //
    //! clear BD count and caches associated with this BD, if no other is using this BD.
    void clearThisBD (bool no_assertion=false);
    //
    //! add the current BD into BD count database
    void addThisBD ();
    //
    struct FabArrayStats
    {
        int  num_fabarrays;
        int  max_num_fabarrays;
        int  max_num_boxarrays;
        int  max_num_ba_use;
        Long num_build;
        FabArrayStats () noexcept : num_fabarrays(0), max_num_fabarrays(0), max_num_boxarrays(0),
                                    max_num_ba_use(1), num_build(0) {;}
        void recordBuild () noexcept {
            ++num_fabarrays;
            ++num_build;
            max_num_fabarrays = std::max(max_num_fabarrays, num_fabarrays);
        }
        void recordDelete () noexcept {
            --num_fabarrays;
        }
        void recordMaxNumBoxArrays (int n) noexcept {
            max_num_boxarrays = std::max(max_num_boxarrays, n);
        }
        void recordMaxNumBAUse (int n) noexcept {
            max_num_ba_use = std::max(max_num_ba_use, n);
        }
        void print () {
            amrex::Print(Print::AllProcs) << "### FabArray ###\n"
                                          << "    tot # of builds       : " << num_build         << "\n"
                                          << "    max # of FabArrays    : " << max_num_fabarrays << "\n"
                                          << "    max # of BoxArrays    : " << max_num_boxarrays << "\n"
                                          << "    max # of BoxArray uses: " << max_num_ba_use    << "\n";
        }
    };
    static FabArrayStats m_FA_stats;

#ifdef BL_USE_MPI
    static bool CheckRcvStats(Vector<MPI_Status>& recv_stats,
			      const Vector<std::size_t>& recv_size,
                              int tag);
#endif

};

std::ostream& operator<< (std::ostream& os, const FabArrayBase::BDKey& id);

Arena* The_FA_Arena ();

}

#endif
