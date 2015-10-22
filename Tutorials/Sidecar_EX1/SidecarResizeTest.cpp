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

	switch(sidecarSignal) {
	  case MySignal:
	    ParallelDescriptor::Bcast(&time_step, 1, 0, ParallelDescriptor::CommunicatorInter());
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << myProcAll << ":: sidecar recv time_step = " << time_step << std::endl;
	    }
	  break;

	  case ParallelDescriptor::SidecarQuitSignal:
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

    BoxLib::Initialize(argc,argv);

    // A flag you need for broadcasting across MPI groups. We always broadcast
    // the data to the sidecar group from the IOProcessor on the compute group.
    int MPI_IntraGroup_Broadcast_Rank;
    int myProcAll(ParallelDescriptor::MyProcAll());
    int nSidecarProcs(0), sidecarSignal(MySignal);
    ParmParse pp;

    pp.get("nSidecars", nSidecarProcs);

    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "nSidecarProcs = " << nSidecarProcs << std::endl;
    }
    ParallelDescriptor::SetNProcsSidecar(nSidecarProcs);

    MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;

    if(ParallelDescriptor::InSidecarGroup()) {
      std::cout << myProcAll << ":: Calling SidecarEventLoop()." << std::endl;
      SidecarEventLoop();
    } else {
      for(int i(0); i < 10; ++i) {  // ----- do time steps.
        sidecarSignal = MySignal;
        ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
        ParallelDescriptor::Bcast(&i, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
      }

      // ---- stop the sidecars
      sidecarSignal = ParallelDescriptor::SidecarQuitSignal;
      ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
    }

    std::cout << myProcAll << ":: Finished first timesteps." << std::endl;
    USleep(1.4);
    ParallelDescriptor::Barrier();

    nSidecarProcs = 4;
    ParallelDescriptor::SetNProcsSidecar(nSidecarProcs);


    MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;

    std::cout << myProcAll << ":: After resizing sidecars:  nProcs = " << ParallelDescriptor::NProcs() << std::endl;

    if(ParallelDescriptor::InSidecarGroup()) {
      std::cout << myProcAll << ":: Second time calling SidecarEventLoop()." << std::endl;
      SidecarEventLoop();
    } else {
      std::cout << myProcAll << ":: Second time doing timesteps." << std::endl;
      for(int i(10); i < 20; ++i) {  // ----- do time steps.
        sidecarSignal = MySignal;
        ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
        ParallelDescriptor::Bcast(&i, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
      }
      USleep(0.4);
      std::cout << myProcAll << ":: Second time after timesteps." << std::endl;

      // ---- stop the sidecars
      sidecarSignal = ParallelDescriptor::SidecarQuitSignal;
      ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
    }


    USleep(0.7);
    ParallelDescriptor::Barrier();

    nSidecarProcs = 0;
    ParallelDescriptor::SetNProcsSidecar(nSidecarProcs);

    USleep(0.3);
    ParallelDescriptor::Barrier();

    std::cout << "_calling Finalize()" << std::endl;
    BoxLib::Finalize();
    return 0;
}
