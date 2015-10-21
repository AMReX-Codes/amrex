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

#include <BoxLib.H>
#include <ParallelDescriptor.H>
#include <Utility.H>
#include <ParmParse.H>

namespace
{
  const unsigned int msps(1000000);
  const int MySignal(42);

  void USleep(double sleepsec) {
    usleep(sleepsec * msps);
  }

  void SidecarEventLoop() {
    bool finished(false);
    int sidecarSignal(-1), time_step(-2);
    int myProcAll(ParallelDescriptor::MyProcAll());

    while ( ! finished) {
        // Receive the signal from the compute group.
        ParallelDescriptor::Bcast(&sidecarSignal, 1, 0, ParallelDescriptor::CommunicatorInter());
        std::cout << myProcAll << ":: sidecar recv signal = " << sidecarSignal << std::endl;

	switch(sidecarSignal) {
	  case MySignal:
	    ParallelDescriptor::Bcast(&time_step, 1, 0, ParallelDescriptor::CommunicatorInter());
            std::cout << myProcAll << ":: sidecar recv time_step = " << time_step << std::endl;
	  break;

	  case ParallelDescriptor::SidecarQuitSignal:
            if (ParallelDescriptor::IOProcessor())
                std::cout << "Sidecars received the quit signal." << std::endl;
            finished = true;
	  break;

	  default:
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

    BoxLib::Initialize(argc,argv);

    int myProcAll(ParallelDescriptor::MyProcAll());
    int nSidecarProcs(0);
    ParmParse pp;

    pp.get("nSidecars", nSidecarProcs);

std::cout << myProcAll << "::_here 0:  nSidecarProcs = " << nSidecarProcs << std::endl;
    ParallelDescriptor::SetNProcsSidecar(nSidecarProcs);

    // A flag you need for broadcasting across MPI groups. We always broadcast
    // the data to the sidecar group from the IOProcessor on the compute group.
    const int MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;

    if(ParallelDescriptor::InSidecarGroup()) {
      SidecarEventLoop();
    } else {

    int sidecarSignal;

    sidecarSignal = MySignal;

//std::cout << myProcAll << "::_here 4." << std::endl;
    // Pretend we're looping over time steps.
    for (unsigned int i = 0; i < 10; ++i)
    {
        ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
        ParallelDescriptor::Bcast(&i, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
    }

    // Don't forget to tell the sidecars to quit! Otherwise they'll keep
    // waiting for signals for more work to do.
    sidecarSignal = ParallelDescriptor::SidecarQuitSignal;
    ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());

    }

    std::cout << "_calling Finalize()" << std::endl;
    BoxLib::Finalize();
    return 0;
}
