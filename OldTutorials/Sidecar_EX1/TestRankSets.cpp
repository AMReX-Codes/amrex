// --------------------------------------------------------------------------
//   TestRankSets.cpp
// --------------------------------------------------------------------------
// An example to test resizing with a disallowed rank set.
// This test is supposed to fail.
// --------------------------------------------------------------------------
#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

using namespace amrex;

// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {

    amrex::Initialize(argc,argv);

    int myProcAll(ParallelDescriptor::MyProcAll());
    int nProcs(ParallelDescriptor::NProcs());
    int nSidecars(1), totalSidecarProcs(4);
    bool printRanks(true);

    if(nProcs < 8) {
      amrex::Abort("**** Error:  must use at least 8 processes.");
    }

    Vector<int> compProcsInAll;
    compProcsInAll.push_back(0);
    for(int i(1); i < nProcs - totalSidecarProcs; ++i) {
      compProcsInAll.push_back(nProcs - totalSidecarProcs - i);  // ---- backwards
    }

    Vector<Vector<int> > sidecarProcsInAll(nSidecars);
    for(int i(nProcs - totalSidecarProcs); i < nProcs; ++i) {
      sidecarProcsInAll[0].push_back(i);
    }
    if(myProcAll == 0) {
      std::cout << myProcAll << ":: Set initial sidecar sizes." << std::endl;
      std::cout << std::endl;
    }
    amrex::USleep(myProcAll);
    ParallelDescriptor::SetNProcsSidecars(compProcsInAll, sidecarProcsInAll, printRanks);


    // ---- now move some procs to the sidecar
    sidecarProcsInAll[0].push_back(compProcsInAll[compProcsInAll.size() - 1]);
    compProcsInAll.pop_back();

    if(myProcAll == 0) {
      std::cout << myProcAll << ":: Move elements from comp to sidecar." << std::endl;
      std::cout << std::endl;
    }
    amrex::USleep(myProcAll);
    ParallelDescriptor::SetNProcsSidecars(compProcsInAll, sidecarProcsInAll, printRanks);


    // ---- now swap the last elements, this test should fail
    int cLast(compProcsInAll[compProcsInAll.size() - 1]);
    compProcsInAll[compProcsInAll.size() - 1] = sidecarProcsInAll[0][sidecarProcsInAll.size() - 1];
    sidecarProcsInAll[0][sidecarProcsInAll.size() - 1] = cLast;

    if(myProcAll == 0) {
      std::cout << myProcAll << ":: Move elements both from and to comp.  This test should fail." << std::endl;
      std::cout << std::endl;
    }
    amrex::USleep(myProcAll);
    ParallelDescriptor::SetNProcsSidecars(compProcsInAll, sidecarProcsInAll, printRanks);

    
    nSidecars = 0;
    ParallelDescriptor::SetNProcsSidecars(nSidecars);

    if(myProcAll == 0) {
      std::cout << myProcAll << ":: Finished." << std::endl;
    }

    amrex::Finalize();
    return 0;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
