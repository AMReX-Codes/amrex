//BL_COPYRIGHT_NOTICE
// -----------------------------------------------------------------
// DistributionMapping.H
// -----------------------------------------------------------------

#ifndef _DISTRIBUTIONMAPPING_H_
#define _DISTRIBUTIONMAPPING_H_

#include <BoxLib.H>
#include <Array.H>
#include <Box.H>
#include <BoxArray.H>

class DistributionMapping {

  public:
    enum DistributionStrategy { ROUNDROBIN, RANDOM, SIZEBALANCED };

  public:
    DistributionMapping();
    DistributionMapping(int nprocessors, const BoxArray &boxes,
		     DistributionStrategy distributionstrategy = ROUNDROBIN);
    ~DistributionMapping();
    void define(int nprocessors, const BoxArray &boxes,
	        DistributionStrategy distributionstrategy = ROUNDROBIN);
    const Array<int> &ProcessorMap() const { return processorMap; }
    int operator()(int level, const Box &box) const;
    int operator[](int index) const { return processorMap[index]; }
    friend ostream& operator<<(ostream &os, const DistributionMapping &pmap);



  private:
    int nProcessors;
    BoxArray boxarray;
    DistributionStrategy distributionStrategy;
    Array<int> processorMap;

    // diagnostics
    Array<int>  objectsPerProcessor;
    Array<long> nPtsPerProcessor;

    void CreateProcessorMap();
};


#endif /*_DISTRIBUTIONMAPPING_H_*/
