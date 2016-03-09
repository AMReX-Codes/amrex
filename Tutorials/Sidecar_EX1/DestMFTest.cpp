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
#include <MultiFab.H>

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
    BoxArray bac, bab;
    MultiFab mf;
    int r(-1);

    while ( ! finished) {
        // Receive the signal from the compute group.
        ParallelDescriptor::Bcast(&sidecarSignal, 1, 0, ParallelDescriptor::CommunicatorInter());

	switch(sidecarSignal) {
	  case MySignal:
        {
	    ParallelDescriptor::Bcast(&time_step, 1, 0, ParallelDescriptor::CommunicatorInter());
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << myProcAll << ":: sidecar recv time_step = " << time_step << std::endl;
	    }
	    bac.clear();
	    BoxArray::RecvBoxArray(bac);
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << myProcAll << ":: sidecar recv ba.size = " << bac.size() << std::endl;
	    }

	    ParallelDescriptor::Bcast(&r, 1, 0, ParallelDescriptor::CommunicatorInter());
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << myProcAll << ":: sidecar recv r = " << r << std::endl;
	    }

        MPI_Group group_sidecar, group_world;
        MPI_Comm_group(ParallelDescriptor::CommunicatorSidecar(), &group_sidecar);
        MPI_Comm_group(MPI_COMM_WORLD, &group_world);

        // Create DM on sidecars with default strategy
        const DistributionMapping dm_sidecar(bac, ParallelDescriptor::NProcsSidecar());
        const Array<int> pm_sidecar = dm_sidecar.ProcessorMap();

        Array<int> pm_world = DistributionMapping::TranslateProcMap(pm_sidecar, group_world, group_sidecar);
        // Don't forget to set the sentinel to the proc # in the new group!
        pm_world[pm_world.size()-1] = ParallelDescriptor::MyProcAll();

        DistributionMapping dm_world(pm_world);
        if (ParallelDescriptor::IOProcessor()) {
          usleep(1000000);
          std::cout << "SIDECAR DM = " << dm_sidecar << std::endl << std::flush;
          std::cout << "WORLD DM = " << dm_world << std::endl << std::flush;
        }

	    //bab = BoxArray(bac).maxSize(r);

	    //mf.clear();
	    //MultiFab::SendMultiFabToSidecars(&mf);
            //if(ParallelDescriptor::IOProcessor()) {
            //  std::cout << myProcAll << ":: sidecar recv mf.ba = " << mf.boxArray() << std::endl;
	    //}

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
    int maxGrid(32), nComp(6), nGhost(2), maxSize(16);
    int ts(0), nSteps(5);
    ParmParse pp;

    pp.query("nSidecars", nSidecarProcs);
    pp.query("maxGrid", maxGrid);
    pp.query("nComp", nComp);
    pp.query("nGhost", nGhost);
    pp.query("maxSize", maxSize);

    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "nSidecarProcs = " << nSidecarProcs << std::endl;
    }

    if(ParallelDescriptor::IOProcessor()) {
      std::cout << myProcAll << ":: Resizing sidecars = " << nSidecarProcs << std::endl;
    }
    ParallelDescriptor::SetNProcsSidecar(nSidecarProcs);
    MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;

    if(ParallelDescriptor::InSidecarGroup()) {
      std::cout << myProcAll << ":: Calling SidecarEventLoop()." << std::endl;
      SidecarEventLoop();
    } else {
      for(int i(ts); i < ts + nSteps; ++i) {  // ----- do time steps.
	if(ParallelDescriptor::IOProcessor()) {
	  std::cout << myProcAll << ":: Doing timestep = " << i << std::endl;
	}

	// Build a Box to span the problem domain, and split it into a BoxArray
	Box baseBox(IntVect(0,0,0), IntVect(maxGrid-1, maxGrid-1, maxGrid-1));
	BoxArray ba(baseBox);
	ba.maxSize(maxSize);

	// This is the DM for the compute processes.
	DistributionMapping comp_DM;
	comp_DM.define(ba, ParallelDescriptor::NProcsComp());

	MultiFab mf(ba, nComp, nGhost, comp_DM, Fab_allocate);

	if(nSidecarProcs > 0) {
	  sidecarSignal = MySignal;
	  ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
	  ParallelDescriptor::Bcast(&i, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
	  BoxArray::SendBoxArray(ba);

	  int r = ( i%4 + 1) * 4;
	  ParallelDescriptor::Bcast(&r, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());

	  //MultiFab::SendMultiFabToSidecars(&mf);

	}
      }
      ts += nSteps;

      if(nSidecarProcs > 0) {
	// ---- stop the sidecars
	sidecarSignal = ParallelDescriptor::SidecarQuitSignal;
	ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
      }
    }
    
    std::cout << myProcAll << ":: Finished timesteps" << std::endl;
    USleep(1.4);
    ParallelDescriptor::Barrier();

    nSidecarProcs = 0;
    ParallelDescriptor::SetNProcsSidecar(nSidecarProcs);

    USleep(0.3);
    ParallelDescriptor::Barrier();

    std::cout << "_calling Finalize()" << std::endl;
    BoxLib::Finalize();
    return 0;
}
