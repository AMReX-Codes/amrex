// --------------------------------------------------------------------------
//   TestRankSets.cpp
// --------------------------------------------------------------------------
// An example to test resizing with a disallowed rank set.
// This test is supposed to fail.
// --------------------------------------------------------------------------
#include <BoxLib.H>
#include <ParallelDescriptor.H>
#include <Utility.H>


// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {

    BoxLib::Initialize(argc,argv);

    int myProcAll(ParallelDescriptor::MyProcAll());
    int nProcs(ParallelDescriptor::NProcs());
    int ioProcNum(ParallelDescriptor::IOProcessorNumber());
    int nSidecars(0);

    nSidecars = 3;

    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "nSidecars = " << nSidecars << std::endl;
    }

    if(ParallelDescriptor::IOProcessor()) {
      std::cout << myProcAll << ":: Resizing sidecars = " << nSidecars << std::endl;
    }

    Array<int> randomRanks;
    if(ParallelDescriptor::IOProcessor()) {
      BoxLib::UniqueRandomSubset(randomRanks, nProcs, nProcs);
      for(int i(0); i < randomRanks.size(); ++i) {
        if(randomRanks[i] == 0) {  // ---- comprank[0] must be 0
	  randomRanks[i] = randomRanks[0];
	  randomRanks[0] = 0;
	}
      }
    }
    BoxLib::BroadcastArray(randomRanks, myProcAll, ioProcNum, ParallelDescriptor::Communicator());

    int totalSidecarProcs(6);
    Array<int> compProcsInAll;

    for(int i(0); i < nProcs - totalSidecarProcs; ++i) {
      compProcsInAll.push_back(randomRanks[i]);
    }

    int sCount(nProcs - totalSidecarProcs);
    Array<Array<int> > sidecarProcsInAll(nSidecars);
    sidecarProcsInAll[0].push_back(randomRanks[sCount++]);
    sidecarProcsInAll[0].push_back(randomRanks[sCount++]);
    sidecarProcsInAll[0].push_back(randomRanks[sCount++]);

    sidecarProcsInAll[1].push_back(randomRanks[sCount++]);

    sidecarProcsInAll[2].push_back(randomRanks[sCount++]);
    sidecarProcsInAll[2].push_back(randomRanks[sCount++]);



    ParallelDescriptor::SetNProcsSidecars(compProcsInAll, sidecarProcsInAll);

    
    ParallelDescriptor::Barrier();
    BoxLib::USleep(myProcAll);
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << myProcAll << ":: Finished timesteps" << std::endl;
    }

    ParallelDescriptor::Barrier();
    nSidecars = 0;
    ParallelDescriptor::SetNProcsSidecars(nSidecars);


    BoxLib::Finalize();
    return 0;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
