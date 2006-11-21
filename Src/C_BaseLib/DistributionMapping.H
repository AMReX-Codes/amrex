#ifndef BL_DISTRIBUTIONMAPPING_H
#define BL_DISTRIBUTIONMAPPING_H

#include <vector>

#include <Pointers.H>
#include <BoxLib.H>
#include <Array.H>
#include <Box.H>

class BoxArray;

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
    //@ManDoc: The distribution strategy: ROUNDROBIN, KNAPSACK or SFC.
    //
    enum Strategy { ROUNDROBIN, KNAPSACK, SFC };
    //
    //@ManDoc: The default constructor.
    //
    DistributionMapping ();
    //
    //@ManDoc: Create object with specified mapping.  The sentinel is already included.
    //
    DistributionMapping (const Array<int>& pmap);
    //
    //@ManDoc: Build mapping out of BoxArray over nprocs processors.
    //
    DistributionMapping (const BoxArray& boxes, int nprocs);
    //
    // This is a very specialized distribution map.
    // Do NOT use it unless you really understand what it does.
    //
    DistributionMapping (const DistributionMapping& d1,
                         const DistributionMapping& d2);
    //
    //@ManDoc: The destructor.
    //
    ~DistributionMapping ();

    /*@ManDoc: Build mapping out of BoxArray over nprocs processors.
               You need to call this if you built your DistributionMapping
               with the default constructor.
    */
    void define (const BoxArray& boxes, int nprocs);

    /*@ManDoc: Returns a constant reference to the mapping of boxes in the
               underlying BoxArray to the CPU that holds the FAB on that Box.
               ProcessorMap()[i] is an integer in the interval [0, NCPU) where
               NCPU is the number of CPUs being used.
    */
    const Array<int>& ProcessorMap () const;
    //
    // Length of the underlying processor map.
    //
    int size () const;
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
               before new MultiFabs are alloc()d.
    */
    static void FlushCache ();
    //
    //@ManDoc: The size of the cache.
    //
    static int CacheSize ();
    //
    //@ManDoc: Append the ProcMap to the Cache.  Checks for consistency.
    //
    static void AddToCache (const DistributionMapping& dm);
    //
    //@ManDoc: Output some simple cache statistics.
    //
    static void CacheStats (std::ostream& os);
    //
    //@ManDoc: Are the distributions equal?
    //
    bool operator== (const DistributionMapping& rhs) const;
    //
    //@ManDoc: Are the distributions different?
    //
    bool operator!= (const DistributionMapping& rhs) const;

    void SFCProcessorMap        (const BoxArray& boxes, const std::vector<long>& wgts, int nprocs);
    void KnapSackProcessorMap   (const std::vector<long>& wgts, int nprocs);
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

protected:
    //
    // Ways to create the processor map.
    //
    void RoundRobinProcessorMap (const BoxArray& boxes, int nprocs);
    void KnapSackProcessorMap   (const BoxArray& boxes, int nprocs);
    void SFCProcessorMap        (const BoxArray& boxes, int nprocs);
    //
    // One which CPU to begin a distribution?
    //
    static int WhereToStart (int nprocs);
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
    // Corresponds to the one specified by `m_Strategy'.
    //
    static PVMF m_BuildMap;
    //
    // Have we initialized from ParmParse yet?
    //
    static bool m_Initialized;

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
        // The length is always equal to boxes.length()+1 where `boxes' is
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
    static std::vector< LnClassPtr<Ref> > m_Cache;
};

inline
int
DistributionMapping::size () const
{
    return m_ref->m_pmap.size();
}

inline
int
DistributionMapping::operator[] (int index) const
{
    return m_ref->m_pmap[index];
}

//
// Our output operator.
//
std::ostream& operator<< (std::ostream& os, const DistributionMapping& pmap);

#endif /*BL_DISTRIBUTIONMAPPING_H*/
