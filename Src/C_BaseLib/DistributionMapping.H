#ifndef BL_DISTRIBUTIONMAPPING_H
#define BL_DISTRIBUTIONMAPPING_H

#include <map>

#include <Pointers.H>
#include <BoxLib.H>
#include <Array.H>
#include <Box.H>

class BoxArray;

//
// Calculates the distribution of FABs to processors.
//
//  This class calculates the distribution of FABs to processors in a
//  FabArray in a multi-processor environment.  By distribution is meant what
//  CPU in the multi-processor environment owns what FAB.  Only the BoxArray
//  on which the FabArray is built is used in determining the distribution.
//  The two types of distributions supported are round-robin and knapsack.  In
//  the round-robin distribution FAB i is owned by CPU i%N where N is total
//  number of CPUs.  In the knapsack distribution the FABs are partitioned
//  across CPUs such that the total volume of the Boxes in the underlying
//  BoxArray are as equal across CPUs as is possible.
//

class DistributionMapping
{
  public:
    //
    // The distribution strategy: ROUNDROBIN, KNAPSACK or SFC.
    //
    enum Strategy { ROUNDROBIN, KNAPSACK, SFC };
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
    DistributionMapping (const Array<int>& pmap, bool put_in_cache = false);
    //
    // Build mapping out of BoxArray over nprocs processors.
    //
    DistributionMapping (const BoxArray& boxes, int nprocs);
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
    void define (const BoxArray& boxes, int nprocs);

    //
    // Build mapping given a pmap.
    //
    void define (const Array<int>& pmap, bool put_in_cache);

    //
    // Returns a constant reference to the mapping of boxes in the
    // underlying BoxArray to the CPU that holds the FAB on that Box.
    // ProcessorMap()[i] is an integer in the interval [0, NCPU) where
    // NCPU is the number of CPUs being used.
    //
    const Array<int>& ProcessorMap () const;
    //
    // Length of the underlying processor map.
    //
    int size () const { return m_ref->m_pmap.size(); }
    //
    // Equivalent to ProcessorMap()[index].
    //
    int operator[] (int index) const { return m_ref->m_pmap[index]; }
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
    // before new MultiFabs are alloc()d.  This removed unreferenced
    // DistributionMaps from the Cache.
    //
    static void FlushCache ();
    //
    // The size of the cache.
    //
    static int CacheSize ();
    //
    // Output some simple cache statistics.
    //
    static void CacheStats (std::ostream& os);
    //
    // Are the distributions equal?
    //
    bool operator== (const DistributionMapping& rhs) const;
    //
    // Are the distributions different?
    //
    bool operator!= (const DistributionMapping& rhs) const;

    void SFCProcessorMap        (const BoxArray& boxes, const std::vector<long>& wgts, int nprocs);
    void KnapSackProcessorMap   (const std::vector<long>& wgts, int nprocs, double* efficiency = 0, bool do_full_knapsack = true);
    void RoundRobinProcessorMap (int nboxes, int nprocs);
    //
    // Initializes distribution strategy from ParmParse.
    //
    // ParmParse options are:
    //
    //   DistributionMapping.strategy = ROUNDROBIN
    //   DistributionMapping.strategy = KNAPSACK
    //   DistributionMapping.strategy = SFC
    //
    static void Initialize ();

    static void Finalize ();

    static bool SameRefs (const DistributionMapping& lhs,
                          const DistributionMapping& rhs) { return lhs.m_ref == rhs.m_ref; }

private:
    //
    // Ways to create the processor map.
    //
    void RoundRobinProcessorMap (const BoxArray& boxes, int nprocs);
    void KnapSackProcessorMap   (const BoxArray& boxes, int nprocs);
    void SFCProcessorMap        (const BoxArray& boxes, int nprocs);

    typedef std::pair<long,int> LIpair;

    struct LIpairComp
    {
        bool operator () (const LIpair& lhs,
                          const LIpair& rhs) const
            {
                return lhs.first < rhs.first;
            }
    };

    static void Sort (std::vector<LIpair>& vec, bool reverse);

    void RoundRobinDoIt (int                  nboxes,
                         int                  nprocs,
                         std::vector<LIpair>* LIpairV = 0);

    void KnapSackDoIt (const std::vector<long>& wgts,
                       int                      nprocs,
                       double&                  efficiency,
                       bool                     do_full_knapsack);

    void SFCProcessorMapDoIt (const BoxArray&          boxes,
                              const std::vector<long>& wgts,
                              int                      nprocs);
    //
    // Least used ordering of CPUs (by # of bytes of FAB data).
    //
    static void LeastUsedCPUs (int nprocs, Array<int>& result);
    //
    // Look for a cached processor map.
    //
    bool GetMap (const BoxArray& boxes);
    //
    // A useful typedef.
    //
    typedef void (DistributionMapping::*PVMF)(const BoxArray&,int);
    //
    // Everyone uses the same Strategy -- defaults to KNAPSACK.
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
    //
    // Our cache of processor maps.
    //
    static std::map< int,LnClassPtr<Ref> > m_Cache;
};

//
// Our output operator.
//
std::ostream& operator<< (std::ostream& os, const DistributionMapping& pmap);

#endif /*BL_DISTRIBUTIONMAPPING_H*/
