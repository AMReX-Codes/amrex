//BL_COPYRIGHT_NOTICE

#ifndef BL_DISTRIBUTIONMAPPING_H
#define BL_DISTRIBUTIONMAPPING_H

//
// $Id: DistributionMapping.H,v 1.11 1998-07-08 16:36:12 lijewski Exp $
//
#ifdef BL_USE_NEW_HFILES
#include <vector>
using std::vector;
#else
#include <vector.h>
#endif

#include <BoxLib.H>
#include <Array.H>
#include <Box.H>
#include <BoxArray.H>

//
//@Man:
//@Memo: Calculates the distribution of FABs to processors.
/*@Doc:

  This class calculates the distribution of FABs to processors in a
  FabArray in a multi-processor environment.  By distribution is meant what
  CPU in the multi-processor environment owns what FAB.  Only the BoxArray
  on which the FabArray is built is used in determining the distribution.
  The two types of distributions supported are round-robin and knapsack.  In
  the round-robin distribution FAB `i' is owned by CPU `i%N' where N is total
  number of CPUs.  In the knapsack distribution the FABs are partitioned
  across CPUs such that the total volume of the Boxes in the underlying
  BoxArray are as equal across CPUs as is possible.
*/

class DistributionMapping
{
  public:
    //
    //@ManDoc: The distribution strategy: ROUNDROBIN or KNAPSACK.
    //
    enum Strategy { ROUNDROBIN, KNAPSACK };
    //
    //@ManDoc: The default constructor.
    //
    DistributionMapping ();
    //
    //@ManDoc: Build mapping on nprocs out of BoxArray.
    //
    DistributionMapping (int             nprocs,
                         const BoxArray& boxes);
    //
    //@ManDoc: The destructor.
    //
    ~DistributionMapping ();

    /*@ManDoc: Build mapping on nprocs out of BoxArray.  You need to
               call this if you built your DistributionMapping with
               the default constructor.
    */
    void define (int             nprocs,
                 const BoxArray& boxes);

    /*@ManDoc: Returns a constant reference to the mapping of boxes in the
               underlying BoxArray to the CPU that holds the FAB on that Box.
               ProcessorMap()[i] is an integer in the interval [0, NCPU) where
               NCPU is the number of CPUs being used.
    */
    const Array<int>& ProcessorMap () const;
    //
    //@ManDoc: Equivalent to ProcessorMap()[index].
    //
    int operator[] (int index) const;
    //
    //@ManDoc: Set distribution strategy.
    //
    static void strategy (Strategy how);
    //
    //@ManDoc: Returns the distribution strategy.
    //
    static Strategy strategy ();

    /*@ManDoc: Flush the cache of processor maps.  The processor map cache
               is only flushed manually.  Only call this after a regridding
               before new MultiFabs are alloced.
    */
    static void FlushCache ();
    //
    //@ManDoc: The size of the caceh.
    //
    static int CacheSize ();
    //
    //@ManDoc: Output some simple cache statistics.
    //
    static void CacheStats (ostream& os);
    //
    //@ManDoc: Are the distributions equal?
    //
    bool operator== (const DistributionMapping& rhs) const;
    //
    //@ManDoc: Are the distributions different?
    //
    bool operator!= (const DistributionMapping& rhs) const;

  private:
    //
    // Two ways to create the processor map.
    //
    void RoundRobinProcessorMap   (int nprocs, const BoxArray& boxes);
    void KnapSackProcessorMap     (int nprocs, const BoxArray& boxes);
    //
    // Initializes distribution strategy from ParmParse.
    //
    // ParmParse options are:
    //
    //   DistributionMapping.strategy = ROUNDROBIN
    //   DistributionMapping.strategy = KNAPSACK
    //
    static void init ();
    //
    // Our cache of processor maps.
    //
    static vector< Array<int> > m_Cache;
    //
    // Look for a cached processor map.
    //
    bool GetMap (int             nprocs,
                 const BoxArray& boxes);
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
    // Local data -- our processor map.
    //
    // The length is always equal to boxes.length()+1 where `boxes' is
    // the BoxArray on which the distribution is based.  It is also true
    // that m_procmap[boxes.length()] == ParallelDescriptor::MyProc().
    // This latter acts as a sentinel in some FabArray loops.
    //
    Array<int> m_procmap;
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
int
DistributionMapping::CacheSize ()
{
    return m_Cache.size();
}

#endif /*BL_DISTRIBUTIONMAPPING_H*/
