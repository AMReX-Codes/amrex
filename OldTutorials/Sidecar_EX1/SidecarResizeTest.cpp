// --------------------------------------------------------------------------
//   SidecarResizeTest.cpp
// --------------------------------------------------------------------------
// An example to test resizing the sidecars and to test the event loop.
// --------------------------------------------------------------------------
#include <new>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sstream>
#include <unistd.h>

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

const int SidecarQuitSignal = -57;

namespace
{
  const unsigned int msps(1000000);
  const int MySignal(42);

//  void USleep(double sleepsec) {
//    usleep(sleepsec * msps);
//  }

  void SidecarEventLoop() {
    bool finished(false);
    int sidecarSignal(-1), time_step(-2);
    int myProcAll(ParallelDescriptor::MyProcAll());

    while ( ! finished) {
        // Receive the signal from the compute group.
        ParallelDescriptor::Bcast(&sidecarSignal, 1, 0, ParallelDescriptor::CommunicatorInter(0));

	switch(sidecarSignal) {
	  case MySignal:
	    ParallelDescriptor::Bcast(&time_step, 1, 0, ParallelDescriptor::CommunicatorInter(0));
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << myProcAll << ":: sidecar recv time_step = " << time_step << std::endl;
	    }
	  break;

	  case SidecarQuitSignal:
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecars received the quit signal." << std::endl;
	    }
            finished = true;
	  break;

	  default:
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "**** Sidecars received bad signal = " << sidecarSignal << std::endl;
	    }
	  break;
	}

    }
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "===== SIDECARS DONE. EXITING ... =====" << std::endl;
    }
  }

}



// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {

    amrex::Initialize(argc,argv);

    // A flag you need for broadcasting across MPI groups. We always broadcast
    // the data to the sidecar group from the IOProcessor on the compute group.
    int MPI_IntraGroup_Broadcast_Rank;
    int myProcAll(ParallelDescriptor::MyProcAll());
    int nSidecarProcs(0), nSidecarProcsFromParmParse(-3), sidecarSignal(MySignal);
    int ts(0), nSteps(10);
    ParmParse pp;

    pp.get("nSidecars", nSidecarProcsFromParmParse);

    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "nSidecarProcs from parmparse = " << nSidecarProcsFromParmParse << std::endl;
    }

    Array<int> howManySidecars(3);
    howManySidecars[0] = nSidecarProcsFromParmParse;
    howManySidecars[1] = 1;
    howManySidecars[2] = 4;

    for(int hMS(0); hMS < howManySidecars.size(); ++hMS) {
      nSidecarProcs = howManySidecars[hMS];
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << myProcAll << ":: Resizing sidecars = " << nSidecarProcs << std::endl;
      }
      ParallelDescriptor::SetNProcsSidecars(nSidecarProcs);
      MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;

      if(ParallelDescriptor::InSidecarGroup()) {
        std::cout << myProcAll << ":: Calling SidecarEventLoop()." << std::endl;
        SidecarEventLoop();
      } else {
        for(int i(ts); i < ts + nSteps; ++i) {  // ----- do time steps.
          if(ParallelDescriptor::IOProcessor()) {
            std::cout << myProcAll << ":: Doing timestep = " << i << std::endl;
          }
	  if(nSidecarProcs > 0) {
            sidecarSignal = MySignal;
            ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter(0));
            ParallelDescriptor::Bcast(&i, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter(0));
	  }
        }
	ts += nSteps;

	if(nSidecarProcs > 0) {
          // ---- stop the sidecars
          sidecarSignal = SidecarQuitSignal;
          ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter(0));
	}
      }

      std::cout << myProcAll << ":: Finished timesteps for hMS = " << hMS << std::endl;
      amrex::USleep(1.4);
      ParallelDescriptor::Barrier();

    }


    nSidecarProcs = 0;
    ParallelDescriptor::SetNProcsSidecars(nSidecarProcs);

    amrex::USleep(0.3);
    ParallelDescriptor::Barrier();

    std::cout << "_calling Finalize()" << std::endl;
    amrex::Finalize();
    return 0;
}
