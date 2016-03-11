// --------------------------------------------------------------------------
//   DestMFTest.cpp
// --------------------------------------------------------------------------
// An example to test copying data from one fabarray to another
//   in another group
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


int nComp(6), nGhost(2);


namespace
{
  const int S_SendBoxArray(42), S_CFATests(43), S_CopyFabArray(44);


  void CopyFabArray(MultiFab *mfSource, MultiFab *mfDest) {
    ParallelDescriptor::Barrier(ParallelDescriptor::CommunicatorAll());

    MultiFab::copyInter(mfSource, mfDest, 0, 0, nComp, 0, 0,
                        ParallelDescriptor::CommunicatorComp(),
                        ParallelDescriptor::CommunicatorSidecar(),
                        ParallelDescriptor::CommunicatorInter(),
                        ParallelDescriptor::InCompGroup());
  }




  void CFATests(BoxArray &ba, DistributionMapping &dm) {

    int myProcAll(ParallelDescriptor::MyProcAll());
    int myProcComp(ParallelDescriptor::MyProcComp());
    int myProcSidecar(ParallelDescriptor::MyProcSidecar());
    MPI_Group group_sidecar(MPI_GROUP_NULL), group_comp(MPI_GROUP_NULL), group_all(MPI_GROUP_NULL);

    BoxLib::USleep(myProcAll / 10.0);
    std::cout << ":::: _in CFATests:  myProcAll myProcComp myProcSidecar = " << myProcAll
              << "  " << myProcComp << "  " << myProcSidecar  << std::endl;
    //std::cout << myProcAll << ":::: _in CFATests:  ba = " << ba << std::endl;
    //std::cout << myProcAll << ":::: _in CFATests:  dm = " << dm << std::endl;

ParallelDescriptor::Barrier(ParallelDescriptor::CommunicatorAll());

    Array<int> pm_sidecar, pm_comp, pm_sidecar_all, pm_comp_all;
    DistributionMapping dm_comp_all, dm_sidecar_all;
    BL_MPI_REQUIRE( MPI_Comm_group(ParallelDescriptor::CommunicatorAll(), &group_all) );

ParallelDescriptor::Barrier(ParallelDescriptor::CommunicatorAll());
    if(ParallelDescriptor::InSidecarGroup()) {
      MPI_Comm_group(ParallelDescriptor::CommunicatorSidecar(), &group_sidecar);
      pm_sidecar = dm.ProcessorMap();
      pm_sidecar_all = DistributionMapping::TranslateProcMap(pm_sidecar, group_all, group_sidecar);
      // Don't forget to set the sentinel to the proc # in the new group!
      pm_sidecar_all[pm_sidecar_all.size()-1] = ParallelDescriptor::MyProcAll();

      MultiFab mfCompTest(ba, 1, 0);
      if(myProcAll == 7) {
        std::cout << myProcSidecar << "::mfCompTest.Dmap = " << mfCompTest.DistributionMap() << std::endl;
      }
    }
ParallelDescriptor::Barrier(ParallelDescriptor::CommunicatorAll());
    if(ParallelDescriptor::InCompGroup()) {
      MPI_Comm_group(ParallelDescriptor::CommunicatorComp(), &group_comp);
      pm_comp = dm.ProcessorMap();
      pm_comp_all = DistributionMapping::TranslateProcMap(pm_comp, group_all, group_comp);
      // Don't forget to set the sentinel to the proc # in the new group!
      pm_comp_all[pm_comp_all.size()-1] = ParallelDescriptor::MyProcAll();

    }

ParallelDescriptor::Barrier(ParallelDescriptor::CommunicatorAll());
    if(ParallelDescriptor::InSidecarGroup()) {
      int MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;
      ParallelDescriptor::Bcast(pm_sidecar_all.dataPtr(), pm_sidecar_all.size(), MPI_IntraGroup_Broadcast_Rank,
	                        ParallelDescriptor::CommunicatorInter());
      pm_comp_all.resize(pm_sidecar_all.size(), -19);  // ---- broadcast these
      ParallelDescriptor::Bcast(pm_comp_all.dataPtr(), pm_comp_all.size(), 0,
	                        ParallelDescriptor::CommunicatorInter());
    }
    if(ParallelDescriptor::InCompGroup()) {
      int MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;
      pm_sidecar_all.resize(pm_comp_all.size(), -17);  // ---- broadcast these
      ParallelDescriptor::Bcast(pm_sidecar_all.dataPtr(), pm_sidecar_all.size(), 0,
	                        ParallelDescriptor::CommunicatorInter());
      ParallelDescriptor::Bcast(pm_comp_all.dataPtr(), pm_comp_all.size(), MPI_IntraGroup_Broadcast_Rank,
	                        ParallelDescriptor::CommunicatorInter());
    }

      dm_sidecar_all.define(pm_sidecar_all, false);
      dm_comp_all.define(pm_comp_all, false);
      BoxLib::USleep(myProcAll / 10.0);
      if(myProcAll == 0 || myProcAll == 7) {
        std::cout << myProcAll << ":::: _in CFATests:  dm_sidecar_all = " << dm_sidecar_all << std::endl;
        std::cout << myProcAll << ":::: _in CFATests:  dm_comp_all = " << dm_comp_all << std::endl;
      }

ParallelDescriptor::Barrier(ParallelDescriptor::CommunicatorAll());
    BoxArray baDest(ba), baSource(ba);
    DistributionMapping dmDest(dm_sidecar_all), dmSource(dm_comp_all);
    MultiFab::CPC cpc(baDest, baSource, dmDest, dmSource);
  }


  void SidecarEventLoop() {
    bool finished(false);
    int sidecarSignal(-1), time_step(-2);
    int myProcAll(ParallelDescriptor::MyProcAll());
    BoxArray bac, bab;
    DistributionMapping sc_DM;
    MultiFab mf;

    while ( ! finished) {  // ---- Receive the signal from the compute group.
        ParallelDescriptor::Bcast(&sidecarSignal, 1, 0, ParallelDescriptor::CommunicatorInter());

	switch(sidecarSignal) {

	  case S_SendBoxArray:
          {
	    bac.clear();
	    BoxArray::RecvBoxArray(bac);
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << myProcAll << ":: sidecar recv ba.size = " << bac.size() << std::endl;
	    }

            MPI_Group group_sidecar, group_all;
            MPI_Comm_group(ParallelDescriptor::CommunicatorSidecar(), &group_sidecar);
            MPI_Comm_group(ParallelDescriptor::CommunicatorAll(), &group_all);

            // Create DM on sidecars with default strategy
            const DistributionMapping dm_sidecar(bac, ParallelDescriptor::NProcsSidecar());
            const Array<int> pm_sidecar = dm_sidecar.ProcessorMap();

            Array<int> pm_all = DistributionMapping::TranslateProcMap(pm_sidecar, group_all, group_sidecar);
            // Don't forget to set the sentinel to the proc # in the new group!
            pm_all[pm_all.size()-1] = ParallelDescriptor::MyProcAll();

            DistributionMapping dm_all(pm_all);
            if (ParallelDescriptor::IOProcessor()) {
              BoxLib::USleep(1);
              std::cout << "SIDECAR DM = " << dm_sidecar << std::endl << std::flush;
              std::cout << "WORLD DM = " << dm_all << std::endl << std::flush;
            }

	    //bab = BoxArray(bac).maxSize(r);

	    //mf.clear();
	    //MultiFab::SendMultiFabToSidecars(&mf);
            //if(ParallelDescriptor::IOProcessor()) {
            //  std::cout << myProcAll << ":: sidecar recv mf.ba = " << mf.boxArray() << std::endl;
	    //}
          }
	  break;

	  case S_CopyFabArray:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecars received the S_CopyFabArray signal." << std::endl;
	    }
	    bac.clear();
	    BoxArray::RecvBoxArray(bac);
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << myProcAll << ":: sidecar recv ba.size = " << bac.size() << std::endl;
	    }
	    sc_DM.define(bac, ParallelDescriptor::NProcsSidecar());
	    mf.define(bac, nComp, nGhost, sc_DM, Fab_allocate);
	    MultiFab *mfSource = 0;
	    MultiFab *mfDest = &mf;
            CopyFabArray(mfSource, mfDest);
          }
	  break;

	  case S_CFATests:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecars received the S_CFATests signal." << std::endl;
	    }
	    sc_DM.define(bac, ParallelDescriptor::NProcsSidecar());
            CFATests(bac, sc_DM);
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
    int nSidecarProcs(0), sidecarSignal(S_SendBoxArray);
    int maxGrid(32), maxSize(16);
    int ts(0), nSteps(5);
    ParmParse pp;

    pp.query("nSidecars", nSidecarProcs);
    pp.query("maxGrid", maxGrid);
    pp.query("nComp", nComp);
    pp.query("nGhost", nGhost);
    pp.query("maxSize", maxSize);
    pp.query("nSteps", nSteps);

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
	  //sidecarSignal = S_SendBoxArray;
	  //ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
	                            //ParallelDescriptor::CommunicatorInter());
	  //BoxArray::SendBoxArray(ba);

	  //MultiFab::SendMultiFabToSidecars(&mf);

	  //sidecarSignal = S_CFATests;
	  sidecarSignal = S_CopyFabArray;
	  ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
	                            ParallelDescriptor::CommunicatorInter());
	  BoxArray::SendBoxArray(ba);
          //CFATests(ba, comp_DM);
	  MultiFab *mfSource = &mf;
	  MultiFab *mfDest = 0;
          CopyFabArray(mfSource, mfDest);

	}
      }
      ts += nSteps;

      if(nSidecarProcs > 0) {
	// ---- stop the sidecars
	sidecarSignal = ParallelDescriptor::SidecarQuitSignal;
	ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
	                          ParallelDescriptor::CommunicatorInter());
      }
    }
    
    ParallelDescriptor::Barrier();
    BoxLib::USleep(myProcAll / 10.0);
    std::cout << myProcAll << ":: Finished timesteps" << std::endl;

    ParallelDescriptor::Barrier();
    nSidecarProcs = 0;
    ParallelDescriptor::SetNProcsSidecar(nSidecarProcs);

    ParallelDescriptor::Barrier();
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "_calling Finalize()" << std::endl;
    }

    BoxLib::Finalize();
    return 0;
}
