//BL_COPYRIGHT_NOTICE
// -----------------------------------------------------------------------------
// DistributionMapping.C
// -----------------------------------------------------------------------------
#include <DistributionMapping.H>
#include <ParallelDescriptor.H>

// -----------------------------------------------------------------------------
DistributionMapping::DistributionMapping()
		 : nProcessors(0),
		   boxarray(),
		   distributionStrategy(ROUNDROBIN),
                   processorMap(),
                   objectsPerProcessor(),
                   nPtsPerProcessor()
{
  CreateProcessorMap();
}

// -----------------------------------------------------------------------------
DistributionMapping::DistributionMapping(int nprocessors, const BoxArray &boxes,
				   DistributionStrategy distributionstrategy)
		 : nProcessors(nprocessors),
		   boxarray(boxes),
		   distributionStrategy(distributionstrategy),
                   processorMap(boxes.length()),
                   objectsPerProcessor(nprocessors),
                   nPtsPerProcessor(nprocessors)
{
  CreateProcessorMap();
}

// -----------------------------------------------------------------------------
DistributionMapping::~DistributionMapping() {
}


// -----------------------------------------------------------------------------
void DistributionMapping::define(int nprocessors, const BoxArray &boxes,
			         DistributionStrategy distributionstrategy)
{
  nProcessors = nprocessors;
  boxarray = boxes;
  distributionStrategy = distributionstrategy;
  processorMap.resize(boxes.length());
  objectsPerProcessor.resize(nprocessors, 0);  // init to zero
  nPtsPerProcessor.resize(nprocessors);
  CreateProcessorMap();
}


// -----------------------------------------------------------------------------
void DistributionMapping::CreateProcessorMap() {
  int i;
  switch(distributionStrategy) {
    case ROUNDROBIN:
      for(i = 0; i < processorMap.length(); i++) {
	processorMap[i] = i % nProcessors;
        ++objectsPerProcessor[processorMap[i]];
      }
    break;
    case RANDOM:
      cerr << "Error in DistributionMapping:  RANDOM not implemented." << endl;
      ParallelDescriptor::Abort("Exiting.");
    break;
    case SIZEBALANCED:
      cerr << "Error in DistributionMapping:  SIZEBALANCED not implemented."
	   << endl;
      ParallelDescriptor::Abort("Exiting.");
    break;
    default:
      cerr << "Error in DistributionMapping:  bad distributionStrategy" << endl;
      ParallelDescriptor::Abort("Exiting.");
  }
}


// -----------------------------------------------------------------------------
int DistributionMapping::operator()(int level, const Box &box) const {
  return -1;
}


// -----------------------------------------------------------------------------
ostream &operator<<(ostream &os, const DistributionMapping &pmap) {
    int i;
    os << "(DistributionMapping" << endl;
    for(i = 0; i < pmap.processorMap.length(); i++) {
      os << "processorMap[" << i << "] = " << pmap.processorMap[i] << endl;
    }
    os << ")" << endl;

    if(os.fail()) {
        BoxLib::Error("operator<<(ostream &, DistributionMapping &) failed");
    }
    return os;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
