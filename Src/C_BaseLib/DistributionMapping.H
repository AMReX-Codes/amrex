//BL_COPYRIGHT_NOTICE

#ifndef BL_DISTRIBUTIONMAPPING_H
#define BL_DISTRIBUTIONMAPPING_H

//
// $Id: DistributionMapping.H,v 1.4 1997-11-25 18:47:05 lijewski Exp $
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
    DistributionMapping (int             nprocessors,
                         const BoxArray& boxes);

    ~DistributionMapping ();

    void define (int             nprocessors,
                 const BoxArray& boxes);

    const Array<int>& ProcessorMap () const { return processorMap; }

    int operator[] (int index) const { return processorMap[index]; }
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
    // Creates the processor map.
    //
    void CreateProcessorMap ();
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
    int        nProcessors;
    BoxArray   boxarray;
    Array<int> processorMap;
    //
    // Everyone uses the same Strategy -- defaults to ROUNDROBIN.
    //
    static Strategy distributionStrategy;
    //
    // Have we initialized from ParmParse yet?
    //
    static bool initialized;
};

//
// Our output operator.
//
ostream&
operator<< (ostream&                   os,
            const DistributionMapping& pmap);

#endif /*BL_DISTRIBUTIONMAPPING_H*/
