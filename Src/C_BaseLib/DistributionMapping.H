//BL_COPYRIGHT_NOTICE

#ifndef BL_DISTRIBUTIONMAPPING_H
#define BL_DISTRIBUTIONMAPPING_H

//
// $Id: DistributionMapping.H,v 1.2 1997-09-15 19:40:05 lijewski Exp $
//

#include <BoxLib.H>
#include <Array.H>
#include <Box.H>
#include <BoxArray.H>

class DistributionMapping
{
  public:

    enum DistributionStrategy { ROUNDROBIN, RANDOM, SIZEBALANCED };

    DistributionMapping ();
    DistributionMapping (int                  nprocessors,
                         const BoxArray&      boxes,
                         DistributionStrategy strategy = ROUNDROBIN);

    ~DistributionMapping ();

    void define (int                  nprocessors,
                 const BoxArray&      boxes,
                 DistributionStrategy strategy = ROUNDROBIN);

    const Array<int>& ProcessorMap () const { return processorMap; }

    int operator() (int        level,
                    const Box& box) const;

    int operator[] (int index) const { return processorMap[index]; }

    friend ostream& operator<< (ostream&                   os,
                                const DistributionMapping& pmap);

  private:

    void CreateProcessorMap ();

    int                  nProcessors;
    BoxArray             boxarray;
    DistributionStrategy distributionStrategy;
    Array<int>           processorMap;
    //
    // diagnostics
    //
    Array<int>           objectsPerProcessor;
    Array<long>          nPtsPerProcessor;
};

#endif /*BL_DISTRIBUTIONMAPPING_H*/
