//BL_COPYRIGHT_NOTICE

//
// $Id: DistributionMapping.cpp,v 1.4 1997-11-18 21:45:28 car Exp $
//

#include <DistributionMapping.H>
#include <ParallelDescriptor.H>

DistributionMapping::DistributionMapping ()
    :
    nProcessors(0),
    boxarray(),
    distributionStrategy(ROUNDROBIN),
    processorMap(),
    objectsPerProcessor(),
    nPtsPerProcessor()
{
    CreateProcessorMap();
}

DistributionMapping::DistributionMapping (int                  nprocessors,
                                          const BoxArray&      boxes,
                                          DistributionStrategy strategy)
    :
    nProcessors(nprocessors),
    boxarray(boxes),
    distributionStrategy(strategy),
    processorMap(boxes.length()),
    objectsPerProcessor(nprocessors),
    nPtsPerProcessor(nprocessors)
{
    CreateProcessorMap();
}

DistributionMapping::~DistributionMapping () {}

void
DistributionMapping::define (int                  nprocessors,
                             const BoxArray&      boxes,
                             DistributionStrategy strategy)
{
    nProcessors = nprocessors;
    boxarray = boxes;
    distributionStrategy = strategy;
    processorMap.resize(boxes.length());
    objectsPerProcessor.resize(nprocessors, 0);
    nPtsPerProcessor.resize(nprocessors);
    CreateProcessorMap();
}

void
DistributionMapping::CreateProcessorMap ()
{
    int i;
    switch (distributionStrategy)
    {
    case ROUNDROBIN:
        for (i = 0; i < processorMap.length(); i++)
        {
            processorMap[i] = i % nProcessors;
            ++objectsPerProcessor[processorMap[i]];
        }
        break;
    case RANDOM:
        BoxLib::Error("RANDOM not implemented");
        break;
    case KNAPSACK:
	BoxLib::Error("KNAPSACK not implement");
	break;
    case SIZEBALANCED:
        BoxLib::Error("SIZEBALANCED not implemented");
        break;
    default:
        BoxLib::Error("Bad DistributionStrategy");
    }
}

int
DistributionMapping::operator () (int        level,
                                  const Box& box) const
{
    return -1;
}

ostream&
operator<< (ostream&                   os,
            const DistributionMapping& pmap)
{
    int i;

    os << "(DistributionMapping" << '\n';
    for (i = 0; i < pmap.processorMap.length(); i++)
    {
        os << "processorMap[" << i << "] = " << pmap.processorMap[i] << '\n';
    }
    os << ')' << '\n';

    if (os.fail())
        BoxLib::Error("operator<<(ostream &, DistributionMapping &) failed");

    return os;
}
