//BL_COPYRIGHT_NOTICE

#ifndef BL_DISTRIBUTIONMAPPING_H
#define BL_DISTRIBUTIONMAPPING_H

//
// $Id: DistributionMapping.H,v 1.6 1997-11-28 00:00:45 lijewski Exp $
//

#include <BoxLib.H>
#include <Array.H>
#include <Box.H>
#include <BoxArray.H>

class DistributionMapping
{
  public:

    enum Strategy { ROUNDROBIN, RANDOM, SIZEBALANCED, KNAPSACK };

    DistributionMapping ();
    DistributionMapping (int             nprocs,
                         const BoxArray& boxes);

    ~DistributionMapping ();

    void define (int             nprocs,
                 const BoxArray& boxes);

    const Array<int>& ProcessorMap () const;

    int operator[] (int index) const;
    //
    // Set distribution strategy.
    //
    static void strategy (Strategy how);
    //
    // Returns the distribution strategy.
    //
    static Strategy strategy ();
    //
    // Flush the cache of processor maps.
    //
    static void FlushCache ();
    //
    // Output some simple cache statistics.
    //
    static void CacheStats (ostream& os);
    //
    // Set the maximum # of processor maps to cache.
    //
    static void MaxCacheSize (int n);
    //
    // The maximum # of processor maps to cache.
    //
    static int MaxCacheSize ();
    //
    // The default maximum cache size.
    //
    enum { MaxSize = 10 };

  private:
    //
    // Four ways to create the processor map.
    //
    void RoundRobinProcessorMap   (int nprocs, const BoxArray& boxes);
    void KnapSackProcessorMap     (int nprocs, const BoxArray& boxes);
    void RandomProcessorMap       (int nprocs, const BoxArray& boxes);
    void SizeBalancedProcessorMap (int nprocs, const BoxArray& boxes);
    //
    // Initializes distribution strategy from ParmParse.
    //
    // ParmParse options are:
    //
    //   DistributionMapping.strategy = ROUNDROBIN
    //   DistributionMapping.strategy = KNAPSACK
    //   DistributionMapping.strategy = RANDOM
    //   DistributionMapping.strategy = SIZEBALANCED
    //
    static void init ();
    //
    // Local data -- our processor map.
    //
    Array<int> m_procmap;
    //
    // Our cache of processor maps.
    //
    static List< Array<int> > m_Cache;
    //
    // Look for a cached processor map.
    //
    static bool GetMap (int             nprocs,
                        const BoxArray& boxes,
                        Array<int>&     procmap);
    //
    // A useful typedef.
    //
    typedef void (DistributionMapping::*PVMF)(int,const BoxArray&);
    //
    // Everyone uses the same Strategy -- defaults to ROUNDROBIN.
    //
    static Strategy m_Strategy;
    //
    // Pointer to one of the CreateProcessorMap() functions.
    // Corresponds to the one specified by `m_Strategy'.
    //
    static PVMF m_BuildMap;
    //
    // Have we initialized from ParmParse yet?
    //
    static bool m_Initialized;
    //
    // Holds the current maximum cache size.
    //
    static int m_MaximumCacheSize;
};

//
// Our output operator.
//
ostream&
operator<< (ostream&                   os,
            const DistributionMapping& pmap);

//
// Inlines.
//

inline
const Array<int>&
DistributionMapping::ProcessorMap () const
{
    return m_procmap;
}

inline
int
DistributionMapping::operator[] (int index) const
{
    return m_procmap[index];
}

inline
DistributionMapping::Strategy
DistributionMapping::strategy ()
{
    return DistributionMapping::m_Strategy;
}

inline
void
DistributionMapping::MaxCacheSize (int n)
{
    assert(n > 0);

    DistributionMapping::m_MaximumCacheSize = n;
}

inline
int
DistributionMapping::MaxCacheSize ()
{
    return DistributionMapping::m_MaximumCacheSize;
}

#endif /*BL_DISTRIBUTIONMAPPING_H*/
