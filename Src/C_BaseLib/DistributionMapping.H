//BL_COPYRIGHT_NOTICE

#ifndef BL_DISTRIBUTIONMAPPING_H
#define BL_DISTRIBUTIONMAPPING_H

//
// $Id: DistributionMapping.H,v 1.5 1997-11-25 22:17:51 lijewski Exp $
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

  private:
    //
    // Four ways to create the processor map.
    //
    void RoundRobinProcessorMap ();
    void KnapSackProcessorMap ();
    void RandomProcessorMap ();
    void SizeBalancedProcessorMap ();
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
    // Local data.
    //
    int        m_nprocs;
    BoxArray   m_boxarray;
    Array<int> m_proc_map;
    //
    // A useful typedef.
    //
    typedef void (DistributionMapping::*PVMF)();
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
    return m_proc_map;
}

inline
int
DistributionMapping::operator[] (int index) const
{
    return m_proc_map[index];
}

#endif /*BL_DISTRIBUTIONMAPPING_H*/
