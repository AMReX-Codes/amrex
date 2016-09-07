#ifndef BL_DISTRIBUTIONMAPPING_H
#define BL_DISTRIBUTIONMAPPING_H

#include <map>
#include <limits>

#include <Pointers.H>
#include <BoxLib.H>
#include <Array.H>
#include <Box.H>
#include <REAL.H>
#include <ParallelDescriptor.H>

class BoxArray;
class MultiFab;

//
// Calculates the distribution of FABs to MPI processes.
//
//  This class calculates the distribution of FABs to MPI processes in a
//  FabArray in a multi-processor environment.  By distribution is meant what
//  MPI process in the multi-processor environment owns what FAB.  Only the BoxArray
//  on which the FabArray is built is used in determining the distribution.
//  The three types of distributions supported are round-robin, knapsack, and SFC.
//  In the round-robin distribution FAB i is owned by CPU i%N where N is total
//  number of CPUs.  In the knapsack distribution the FABs are partitioned
//  across CPUs such that the total volume of the Boxes in the underlying
//  BoxArray are as equal across CPUs as is possible.  The SFC distribution is
//  based on a space filling curve.
//

class DistributionMapping
{
  public:
    //
    // The distribution strategies
    //
    enum Strategy { UNDEFINED = -1, ROUNDROBIN, KNAPSACK, SFC, PFC, RRSFC };
    //
    // The default constructor.
    //
    DistributionMapping ();
    //
    // Create an object with the specified mapping.
    // VERY IMPORTANT:  The size of pmap must be (yourboxarray.size() + 1).
    //                  The value in pmap[yourboxarray.size()] is a sentinel used
    //                  in boxlib and its value must be set to
    //                  ParallelDescriptor::MyProc() by the caller.
    //
    DistributionMapping (const Array<int>& pmap, bool put_in_cache = false,
			 ParallelDescriptor::Color color = ParallelDescriptor::DefaultColor());
    //
    // Build mapping out of BoxArray over nprocs processors.
    //
    DistributionMapping (const BoxArray& boxes, int nprocs, 
			 ParallelDescriptor::Color color = ParallelDescriptor::DefaultColor());
    //
    // This is a very specialized distribution map.
    // Do NOT use it unless you really understand what it does.
    //
    DistributionMapping (const DistributionMapping& d1,
                         const DistributionMapping& d2);
    //
    // The destructor.
    //
    ~DistributionMapping ();
    //
    // Copy constructor.
    //
    DistributionMapping (const DistributionMapping& rhs);
    //
    // Assignment operator.
    //
    DistributionMapping& operator= (const DistributionMapping& rhs);
    //
    // Build mapping out of BoxArray over nprocs processors.
    // You need to call this if you built your DistributionMapping
    // with the default constructor.
    //
    void define (const BoxArray& boxes, int nprocs,
		 ParallelDescriptor::Color color = ParallelDescriptor::DefaultColor());
    //
    // Build mapping out of an Array of ints. You need to call this if you
    // built your DistributionMapping with the default constructor.
    //
    void define (const Array<int>& pmap);
    void define (const Array<int>& pmap, bool put_in_cache);
    //
    // Returns a constant reference to the mapping of boxes in the
    // underlying BoxArray to the CPU that holds the FAB on that Box.
    // ProcessorMap()[i] is an integer in the interval [0, NCPU) where
    // NCPU is the number of CPUs being used.
    //
    const Array<int>& ProcessorMap () const;
    //
    ParallelDescriptor::Color color() const { return m_color; }
    //
    // Length of the underlying processor map.
    //
    long size () const { return m_ref->m_pmap.size(); }
    long capacity () const { return m_ref->m_pmap.capacity(); }
    //
    // Number of references to this DistributionMapping
    //
    int linkCount () const { return m_ref.linkCount(); }
    //
    // Equivalent to ProcessorMap()[index].
    //
    int operator[] (int index) const { return m_ref->m_pmap[index]; }
    //
    // Replace the cached processor map.  This is to support FabArray::MoveAllFabs
    //   All FabArrays using the cached map must move their fabs before
    //   calling this function.
    void ReplaceCachedProcessorMap(const Array<int> &newProcMap);
    //
    // Set/get the distribution strategy.
    //
    static void strategy (Strategy how);

    static Strategy strategy ();
    //
    // Set/get the space filling curve threshold.
    //
    static void SFC_Threshold (int n);

    static int SFC_Threshold ();
    //
    // Flush the cache of processor maps.  The processor map cache
    // is only flushed manually.  Only call this after a regridding
    // before new MultiFabs are alloc()d.  This removes unreferenced
    // DistributionMaps from the Cache.
    //
    static void FlushCache ();
    //
    // Delete the cache.  This is to support dynamic sidecars.
    //
    static void DeleteCache ();
    //
    // The size of the cache.
    //
    static int CacheSize ();
    //
    // Output some simple cache statistics.
    //
    static void CacheStats (std::ostream& os, int whichProc = 0);
    //
    // Put in cache
    //
    void PutInCache();
    //
    // Are the distributions equal?
    //
    bool operator== (const DistributionMapping& rhs) const;
    //
    // Are the distributions different?
    //
    bool operator!= (const DistributionMapping& rhs) const;

    void SFCProcessorMap(const BoxArray& boxes, const std::vector<long>& wgts,
                         int nprocs);
    void PFCProcessorMap(const BoxArray& boxes, const std::vector<long>& wgts,
                         int nprocs);
    void KnapSackProcessorMap(const std::vector<long>& wgts, int nprocs,
                              Real* efficiency = 0,
			      bool do_full_knapsack = true,
			      int nmax = std::numeric_limits<int>::max());
    void RoundRobinProcessorMap(int nboxes, int nprocs);
    //
    // Initializes distribution strategy from ParmParse.
    //
    // ParmParse options are:
    //
    //   DistributionMapping.strategy = ROUNDROBIN
    //   DistributionMapping.strategy = KNAPSACK
    //   DistributionMapping.strategy = SFC
    //   DistributionMapping.strategy = PFC
    //   DistributionMapping.strategy = RRFC
    //
    static void Initialize ();

    static void Finalize ();

    static bool SameRefs (const DistributionMapping& lhs,
                          const DistributionMapping& rhs)
		  { return lhs.m_ref == rhs.m_ref; }

    static void PrintDiagnostics(const std::string &filename);
    //
    // Initialize the topological proximity map
    //
    static void InitProximityMap(bool makeMap = false, bool reinit = false);
    static int NHops(const Box &tbox, const IntVect &ivfrom, const IntVect &ivto);
    static int ProcNumberFromRank(const int rank);
    static std::vector<int> RanksFromProcNumber(const int procnum);
    static IntVect TopIVFromProcNumber(const int procnum);
    static std::vector<int> ProcNumbersFromTopIV(const IntVect &iv);
    static IntVect TopIVFromRank(const int rank);
    static std::vector<int> RanksFromTopIV(const IntVect &iv);
    static std::string GetProcName();
    static int GetProcNumber();
    static int ProximityMap(const int rank)   { return proximityMap[rank];   }
    static int ProximityOrder(const int rank) { return proximityOrder[rank]; }

#if !(defined(BL_NO_FORT) || defined(WIN32))
    static void ReadCheckPointHeader(const std::string &filename,
				     Array<IntVect>  &refRatio,
                                     Array<BoxArray> &allBoxes);
#endif
    static void PFCMultiLevelMap(const Array<IntVect>  &refRatio,
                                 const Array<BoxArray> &allBoxes);

    static Array<Array<int> > MultiLevelMapPFC(const Array<IntVect>  &refRatio,
                                               const Array<BoxArray> &allBoxes,
					       int maxgrid);
    static Array<Array<int> > MultiLevelMapRandom(const Array<IntVect>  &refRatio,
                                                  const Array<BoxArray> &allBoxes,
						  int maxgrid,
						  int maxRank = -1, int minRank = 0);
    static Array<Array<int> > MultiLevelMapKnapSack(const Array<IntVect>  &refRatio,
                                                    const Array<BoxArray> &allBoxes,
						    int maxgrid);

    static int NDistMaps() { return nDistMaps; }
    static void SetNDistMaps(int ndm) { nDistMaps = ndm; }
    int DistMapID() const { return dmID; }
    void SetDistMapID(int dmid) { dmID = dmid; }

    bool Check() const;

    //
    // This gives a unique ID of the reference, which is different from dmID above.
    //
    ptrdiff_t getRefID () const;

#ifdef BL_USE_MPI
    static Array<int> TranslateProcMap(const Array<int> &pm_old, const MPI_Group group_new, const MPI_Group group_old);
#endif

    static DistributionMapping makeKnapSack (const MultiFab& weight);

private:
    //
    // Ways to create the processor map.
    //
    void RoundRobinProcessorMap (const BoxArray& boxes, int nprocs);
    void KnapSackProcessorMap   (const BoxArray& boxes, int nprocs);
    void SFCProcessorMap        (const BoxArray& boxes, int nprocs);
    void PFCProcessorMap        (const BoxArray& boxes, int nprocs);
    void RRSFCProcessorMap      (const BoxArray& boxes, int nprocs);

    typedef std::pair<long,int> LIpair;

    struct LIpairLT
    {
        bool operator () (const LIpair& lhs,
                          const LIpair& rhs) const
            {
                return lhs.first < rhs.first;
            }
    };

    struct LIpairGT
    {
        bool operator () (const LIpair& lhs,
                          const LIpair& rhs) const
            {
                return lhs.first > rhs.first;
            }
    };

    static void Sort (std::vector<LIpair>& vec, bool reverse);

    void RoundRobinDoIt (int                  nboxes,
                         int                  nprocs,
                         std::vector<LIpair>* LIpairV = 0);

    void KnapSackDoIt (const std::vector<long>& wgts,
                       int                      nprocs,
                       Real&                    efficiency,
                       bool                     do_full_knapsack,
		       int                      nmax=std::numeric_limits<int>::max());

    void SFCProcessorMapDoIt (const BoxArray&          boxes,
                              const std::vector<long>& wgts,
                              int                      nprocs);

    void PFCProcessorMapDoIt (const BoxArray&          boxes,
                              const std::vector<long>& wgts,
                              int                      nprocs);

    void RRSFCDoIt           (const BoxArray&          boxes,
                              int                      nprocs);

    //
    // Current # of bytes of FAB data.
    //
    static void CurrentBytesUsed (int nprocs, Array<long>& result);
    static void CurrentCellsUsed (int nprocs, Array<long>& result);
    //
    // Least used ordering of CPUs (by # of bytes of FAB data).
    //
    void LeastUsedCPUs (int nprocs, Array<int>& result);
    //
    // rteam: Least used ordering of Teams
    // rworker[i]: Least used ordering of team workers for Team i
    //
    void LeastUsedTeams (Array<int>& rteam, Array<Array<int> >& rworker, int nteams, int nworkers);
    //
    // Look for a cached processor map.
    //
    bool GetMap (const BoxArray& boxes);
    bool GetMap (int nBoxes);
    //
    // A useful typedef.
    //
    typedef void (DistributionMapping::*PVMF)(const BoxArray &, int);
    //
    // Everyone uses the same Strategy -- defaults to SFC.
    //
    static Strategy m_Strategy;
    //
    // Pointer to one of the CreateProcessorMap() functions.
    // Corresponds to the one specified by m_Strategy.
    //
    static PVMF m_BuildMap;

    class Ref
    {
        friend class DistributionMapping;
        //
        // Constructors to match those in DistributionMapping ....
        //
        Ref ();

        Ref (int len);

        Ref (const Array<int>& pmap);

        Ref (const Ref& rhs);
        //
        // Local data -- our processor map.
        //
        // The length is always equal to boxes.length()+1 where boxes is
        // the BoxArray on which the distribution is based.  It is also true
        // that m_pmap[boxes.length()] == ParallelDescriptor::MyProc().
        // This latter acts as a sentinel in some FabArray loops.
        //
        Array<int> m_pmap;
    };
    //
    // The data -- a reference-counted pointer to a Ref.
    //
    LnClassPtr<Ref> m_ref;
    ParallelDescriptor::Color m_color;
    //
    // Our cache of processor maps.
    //
    static std::map< std::pair<int,int>, LnClassPtr<Ref> > m_Cache;
    //
    // Topological proximity map
    //
    static long totalCells;
    static Real bytesPerCell;
    static Array<int> proximityMap;    // i == rank, pMap[i]   == proximity mapped rank
    static Array<int> proximityOrder;  // i == rank, pOrder[i] == proximity mapped order
    static Array<long> totalBoxPoints;  // i == rank

    static int nDistMaps;
    int dmID;
};

//
// Our output operator.
//
std::ostream& operator<< (std::ostream& os, const DistributionMapping& pmap);

#endif /*BL_DISTRIBUTIONMAPPING_H*/
