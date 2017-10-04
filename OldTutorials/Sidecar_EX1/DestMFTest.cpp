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

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>

using namespace amrex;

int nComp(6), nGhost(2);


// --------------------------------------------------------------------------
namespace
{
  const int S_SendBoxArray(42), S_CFATests(43);
  const int S_CopyFabArray(44), S_CopyFabArrayFromSidecar(45);
  const int S_QuitSignal(-100);


  // --------------------------------------------------------------------------
  void CopyFabArray(MultiFab *mfSource, MultiFab *mfDest,
                    int srcComp, int destComp, int numComp,
		    int srcNGhost, int destNGhost,
		    bool fromComp)
  {

    static int count(0);
    std::stringstream css;
    css << "TS_" << count << "_";

    VisMF::SetNOutFiles(1);

    BL_ASSERT( (mfSource == 0) || (mfDest == 0) );
    BL_ASSERT( ! ((mfSource == 0) && (mfDest == 0)) );

    bool isSrc;
    if(mfSource != 0) {
      isSrc = true;
    } else {
      isSrc = false;
    }

    MPI_Comm commSrc, commDest, commInter, commBoth;
    commInter = ParallelDescriptor::CommunicatorInter(0);
    commBoth = ParallelDescriptor::CommunicatorBoth(0);
    if(fromComp) {
      commSrc  = ParallelDescriptor::CommunicatorComp();
      commDest = ParallelDescriptor::CommunicatorSidecar();
    } else {
      commSrc  = ParallelDescriptor::CommunicatorSidecar();
      commDest = ParallelDescriptor::CommunicatorComp();
    }

    if(mfSource != 0) {
      VisMF::Write(*mfSource, css.str() + "mfSource_before");
    }
    if(mfDest != 0) {
      VisMF::Write(*mfDest, css.str() + "mfDest_before");
    }
    MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, numComp,
                        srcNGhost, destNGhost,
                        commSrc, commDest, commInter, commBoth,
                        isSrc);

    if(mfSource != 0) {
      VisMF::Write(*mfSource, css.str() + "mfSource_after");
    }
    if(mfDest != 0) {
      VisMF::Write(*mfDest, css.str() + "mfDest_after");
    }
    ++count;
  }




  // --------------------------------------------------------------------------
  void CFATests(BoxArray &ba, DistributionMapping &dm) {

    int myProcAll(ParallelDescriptor::MyProcAll());
    int myProcComp(ParallelDescriptor::MyProcComp());
    int myProcSidecar(ParallelDescriptor::MyProcSidecar());
    bool addToCache(false);
    MPI_Group group_sidecar(MPI_GROUP_NULL), group_all(MPI_GROUP_NULL);

    ParallelDescriptor::Barrier(ParallelDescriptor::CommunicatorAll());
    amrex::USleep(myProcAll / 10.0);
    std::cout << ":::: _in CFATests:  myProcAll myProcComp myProcSidecar = " << myProcAll
              << "  " << myProcComp << "  " << myProcSidecar  << std::endl;

    Array<int> pm_sidecar, pm_comp, pm_sidecar_all, pm_comp_all;
    DistributionMapping dm_comp_all, dm_sidecar_all;
    BL_MPI_REQUIRE( MPI_Comm_group(ParallelDescriptor::CommunicatorAll(), &group_all) );

    if(ParallelDescriptor::InSidecarGroup()) {
      MPI_Comm_group(ParallelDescriptor::CommunicatorSidecar(), &group_sidecar);
      pm_sidecar = dm.ProcessorMap();
      pm_sidecar_all = DistributionMapping::TranslateProcMap(pm_sidecar, group_all, group_sidecar);

      dm_sidecar_all.define(pm_sidecar_all);
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << myProcAll << ":::: _in CFATests:  dm = " << dm << std::endl;
        std::cout << myProcAll << ":::: _in CFATests:  dm_sidecar_all = " << dm_sidecar_all << std::endl;
      }
    }

      dm_comp_all.define(pm_comp_all);
      amrex::USleep(myProcAll / 10.0);
      if(myProcAll == 0) {
        std::cout << myProcAll << ":::: _in CFATests:  dm_comp = " << dm << std::endl;
      }
  }


  // --------------------------------------------------------------------------
  void SidecarEventLoop() {
    bool finished(false);
    int sidecarSignal(-1);
    int myProcAll(ParallelDescriptor::MyProcAll());
    BoxArray bac, bab;
    DistributionMapping sc_DM;

    while ( ! finished) {  // ---- Receive the signal from the compute group.

        ParallelDescriptor::Bcast(&sidecarSignal, 1, 0, ParallelDescriptor::CommunicatorInter(0));

	switch(sidecarSignal) {

	  case S_SendBoxArray:
          {
	    bac.clear();
	    BoxArray::RecvBoxArray(bac, 0);
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << myProcAll << ":: sidecar recv ba.size = " << bac.size() << std::endl;
	    }

            MPI_Group group_sidecar, group_all;
            MPI_Comm_group(ParallelDescriptor::CommunicatorSidecar(), &group_sidecar);
            MPI_Comm_group(ParallelDescriptor::CommunicatorAll(), &group_all);

            // Create DM on sidecars with default strategy
            const DistributionMapping dm_sidecar(bac, ParallelDescriptor::NProcsSidecar(0));
            const Array<int> pm_sidecar = dm_sidecar.ProcessorMap();

            Array<int> pm_all = DistributionMapping::TranslateProcMap(pm_sidecar, group_all, group_sidecar);

            DistributionMapping dm_all(pm_all);
            if (ParallelDescriptor::IOProcessor()) {
              amrex::USleep(1);
              std::cout << "SIDECAR DM = " << dm_sidecar << std::endl << std::flush;
              std::cout << "WORLD DM = " << dm_all << std::endl << std::flush;
            }
          }
	  break;

	  case S_CopyFabArray:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecars received the S_CopyFabArray signal." << std::endl;
	    }
	    bac.clear();
	    BoxArray::RecvBoxArray(bac, 0);
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << myProcAll << ":: sidecar recv ba.size = " << bac.size() << std::endl;
	    }
	    sc_DM.define(bac, ParallelDescriptor::NProcsSidecar(0));
            MultiFab mf(bac, sc_DM, nComp, nGhost);
	    mf.setVal(-1.0);
	    MultiFab *mfSource = 0;
	    MultiFab *mfDest = &mf;
	    bool fromComp(true);
            CopyFabArray(mfSource, mfDest, 0, 0, nComp, 0, 0, fromComp);
          }
	  break;

	  case S_CopyFabArrayFromSidecar:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecars received the S_CopyFabArrayFromSidecar signal." << std::endl;
	    }
            Box baseBox(IntVect(2,4,6), IntVect(15, 19, 79));
            BoxArray ba(baseBox);
            ba.maxSize(4);
	    DistributionMapping dm{ba};
	    MultiFab mf(ba, dm, nComp, nGhost);
	    for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
	      for(int i(0); i < mf[mfi].nComp(); ++i) {
	        mf[mfi].setVal(myProcAll + (Real) i / 1000.0, i);
	      }
	    }
	    MultiFab *mfSource = &mf;
	    MultiFab *mfDest = 0;
	    int srcComp(4), destComp(2), numComp(1);
	    bool fromComp(false);
            CopyFabArray(mfSource, mfDest, srcComp, destComp, numComp, 0, 0, fromComp);
          }
	  break;

	  case S_CFATests:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecars received the S_CFATests signal." << std::endl;
	    }
	    sc_DM.define(bac, ParallelDescriptor::NProcsSidecar(0));
            CFATests(bac, sc_DM);
          }
	  break;

	  case S_QuitSignal:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecars received the quit signal." << std::endl;
	    }
            finished = true;
	  }
	  break;

	  default:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "**** Sidecars received bad signal = " << sidecarSignal << std::endl;
	    }
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
    ParallelDescriptor::SetNProcsSidecars(nSidecarProcs);
    MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;


    if(ParallelDescriptor::InSidecarGroup()) {

      amrex::USleep(myProcAll / 10.0);
      std::cout << myProcAll << ":: Calling SidecarEventLoop()." << std::endl;
      SidecarEventLoop();

    } else {

      // ---- Build a Box to span the problem domain, and split it into a BoxArray
      Box baseBox(IntVect(0,0,0), IntVect(maxGrid-1, maxGrid-1, maxGrid-1));
      BoxArray ba(baseBox);
      ba.maxSize(maxSize);

      // ---- This is the DM for the compute processes.
      DistributionMapping comp_DM;
      comp_DM.define(ba, ParallelDescriptor::NProcsComp());


      for(int i(ts); i < ts + nSteps; ++i) {  // ----- do time steps.
	if(ParallelDescriptor::IOProcessor()) {
	  std::cout << myProcAll << ":: Doing timestep = " << i << std::endl;
	}

	if(nSidecarProcs > 0) {

	  if((i - ts) == 0) {  // ---- do a simple mf copy test
	    MultiFab mf(ba, comp_DM, nComp, nGhost);
	    for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
	      for(int i(0); i < mf[mfi].nComp(); ++i) {
	        mf[mfi].setVal(myProcAll + (Real) i / 1000.0, i);
	      }
	    }

	    sidecarSignal = S_CopyFabArray;
	    ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
	                              ParallelDescriptor::CommunicatorInter(0));
	    BoxArray::SendBoxArray(ba, 0);

	    MultiFab *mfSource = &mf;
	    MultiFab *mfDest = 0;
	    bool fromComp(true);
            CopyFabArray(mfSource, mfDest, 0, 0, nComp, 0, 0, fromComp);
	  }

	  if((i - ts) == 1) {  // ---- do a shrinked boxarray mf copy test
	    BoxArray bashrink(ba);
	    bashrink.grow(-1);
	    MultiFab mf(bashrink, comp_DM, nComp, nGhost);
	    for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
	      for(int i(0); i < mf[mfi].nComp(); ++i) {
	        mf[mfi].setVal(myProcAll + (Real) i / 1000.0, i);
	      }
	    }

	    sidecarSignal = S_CopyFabArray;
	    ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
	                              ParallelDescriptor::CommunicatorInter(0));
	    BoxArray::SendBoxArray(ba, 0);  // ---- send the sidecar the unshrunk boxarray

	    MultiFab *mfSource = &mf;
	    MultiFab *mfDest = 0;
	    bool fromComp(true);
            CopyFabArray(mfSource, mfDest, 0, 0, nComp, 0, 0, fromComp);
	  }


	  if((i - ts) == 2) {  // ---- do a shifted boxarray mf copy test
	    BoxArray bashift(ba);
	    bashift.shift(IntVect(1,2,3));
	    MultiFab mf(bashift, comp_DM, nComp, nGhost);
	    for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
	      for(int i(0); i < mf[mfi].nComp(); ++i) {
	        mf[mfi].setVal(myProcAll + (Real) i / 1000.0, i);
	      }
	    }

	    sidecarSignal = S_CopyFabArray;
	    ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
	                              ParallelDescriptor::CommunicatorInter(0));
	    BoxArray::SendBoxArray(ba, 0);  // ---- send the sidecar the unshifted boxarray

	    MultiFab *mfSource = &mf;
	    MultiFab *mfDest = 0;
	    bool fromComp(true);
            CopyFabArray(mfSource, mfDest, 0, 0, nComp, 0, 0, fromComp);
	  }


	  if((i - ts) == 3) {  // ---- copy part of a FabArray from the sidecar
	    MultiFab mf(ba, comp_DM, nComp, nGhost);
	    mf.setVal(-2.0);

	    sidecarSignal = S_CopyFabArrayFromSidecar;
	    ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
	                              ParallelDescriptor::CommunicatorInter(0));

	    MultiFab *mfSource = 0;
	    MultiFab *mfDest = &mf;
	    int srcComp(4), destComp(2), numComp(1);
	    bool fromComp(false);
            CopyFabArray(mfSource, mfDest, srcComp, destComp, numComp, 0, 0, fromComp);
	  }


	  // ---- other tests?
	  //sidecarSignal = S_CFATests;
	  //ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
	                            //ParallelDescriptor::CommunicatorInter(0));
          //CFATests(ba, comp_DM);

	}
      }
      ts += nSteps;

      if(nSidecarProcs > 0) {
	// ---- stop the sidecars
	sidecarSignal = S_QuitSignal;
	ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
	                          ParallelDescriptor::CommunicatorInter(0));
      }
    }
    
    ParallelDescriptor::Barrier();
    amrex::USleep(myProcAll / 10.0);
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << myProcAll << ":: Finished timesteps" << std::endl;
    }

    ParallelDescriptor::Barrier();
    nSidecarProcs = 0;
    ParallelDescriptor::SetNProcsSidecars(nSidecarProcs);

    amrex::Finalize();

    return 0;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
