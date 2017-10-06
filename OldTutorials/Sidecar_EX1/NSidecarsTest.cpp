// --------------------------------------------------------------------------
//   NSidecarsTest.cpp
// --------------------------------------------------------------------------
// An example to test multiple sidecars
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
  const int S_SendBoxArray(42);
  const int S_CopyFabArray(44), S_CopyFabArrayFromSidecar(45);
  const int S_ResizeSidecars(46);
  const int S_QuitSignal(-100);


  // --------------------------------------------------------------------------
  void CopyFabArray(MultiFab *mfSource, MultiFab *mfDest,
                    int srcComp, int destComp, int numComp,
		    int srcNGhost, int destNGhost,
		    bool fromComp, int whichSidecar)
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
    commInter = ParallelDescriptor::CommunicatorInter(whichSidecar);
    commBoth  = ParallelDescriptor::CommunicatorBoth(whichSidecar);  // ---- both src and dest
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
  void SidecarEventLoop2() {
    int whichSidecar(2);
    bool finished(false);
    int sidecarSignal(-1);
    int myProcAll(ParallelDescriptor::MyProcAll());
    BoxArray bac, bab;
    DistributionMapping sc_DM;

    while ( ! finished) {  // ---- Receive the signal from the compute group.

        ParallelDescriptor::Bcast(&sidecarSignal, 1, 0,
	                          ParallelDescriptor::CommunicatorInter(whichSidecar));

	switch(sidecarSignal) {

	  case S_SendBoxArray:
          {
	    bac.clear();
	    BoxArray::RecvBoxArray(bac, whichSidecar);
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << myProcAll << ":: sidecar " << whichSidecar << " recv ba.size = "
	                << bac.size() << std::endl;
	    }

            MPI_Group group_sidecar, group_all;
            MPI_Comm_group(ParallelDescriptor::CommunicatorSidecar(), &group_sidecar);
            MPI_Comm_group(ParallelDescriptor::CommunicatorAll(), &group_all);

            // Create DM on sidecars with default strategy
            const DistributionMapping dm_sidecar(bac, ParallelDescriptor::NProcsSidecar(whichSidecar));
            const Vector<int> pm_sidecar = dm_sidecar.ProcessorMap();

            Vector<int> pm_all = DistributionMapping::TranslateProcMap(pm_sidecar, group_all, group_sidecar);

            DistributionMapping dm_all(pm_all);
            if (ParallelDescriptor::IOProcessor()) {
              amrex::USleep(1);
              std::cout << "SIDECAR " << whichSidecar << " DM = " << dm_sidecar << std::endl << std::flush;
              std::cout << "WORLD DM = " << dm_all << std::endl << std::flush;
            }
          }
	  break;

	  case S_CopyFabArray:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecar " << whichSidecar << " received the S_CopyFabArray signal." << std::endl;
	    }
	    bac.clear();
	    BoxArray::RecvBoxArray(bac, whichSidecar);
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << myProcAll << ":: sidecar " << whichSidecar << " recv ba.size = "
	                << bac.size() << std::endl;
	    }
	    sc_DM.define(bac, ParallelDescriptor::NProcsSidecar(whichSidecar));
            MultiFab mf(bac, sc_DM, nComp, nGhost);
	    mf.setVal(-1.0);
	    MultiFab *mfSource = 0;
	    MultiFab *mfDest = &mf;
	    bool fromComp(true);
            CopyFabArray(mfSource, mfDest, 0, 0, nComp, 0, 0, fromComp, whichSidecar);
          }
	  break;

	  case S_CopyFabArrayFromSidecar:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecar " << whichSidecar << " received the S_CopyFabArrayFromSidecar signal."
	                << std::endl;
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
            CopyFabArray(mfSource, mfDest, srcComp, destComp, numComp, 0, 0, fromComp, whichSidecar);
          }
	  break;

	  case S_QuitSignal:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecar " << whichSidecar << " received the quit signal." << std::endl;
	    }
            finished = true;
	  }
	  break;

	  default:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "**** Sidecar " << whichSidecar << " received bad signal = "
	                << sidecarSignal << std::endl;
	    }
	  }
	  break;
	}

    }
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "===== SIDECARS " << whichSidecar << " DONE. EXITING ... =====" << std::endl;
    }
  }


  // --------------------------------------------------------------------------
  void SidecarEventLoop1() {
    int whichSidecar(1);
    bool finished(false);
    int sidecarSignal(-1);

    while ( ! finished) {  // ---- Receive the signal from the compute group.

        ParallelDescriptor::Bcast(&sidecarSignal, 1, 0,
	                          ParallelDescriptor::CommunicatorInter(whichSidecar));

	switch(sidecarSignal) {

	  case S_QuitSignal:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecar 1 received the quit signal." << std::endl;
	    }
            finished = true;
	  }
	  break;

	  default:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "**** Sidecar " << whichSidecar << " received bad signal = "
	                << sidecarSignal << std::endl;
	    }
	  }
	  break;
	}

    }
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "===== SIDECAR " << whichSidecar << " DONE. EXITING ... =====" << std::endl;
    }

  }


  // --------------------------------------------------------------------------
  void SidecarEventLoop0() {
    int whichSidecar(0);
    bool finished(false);
    int sidecarSignal(-1);

    while ( ! finished) {  // ---- Receive the signal from the compute group.

        ParallelDescriptor::Bcast(&sidecarSignal, 1, 0,
	                          ParallelDescriptor::CommunicatorInter(whichSidecar));

	switch(sidecarSignal) {

	  case S_QuitSignal:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "Sidecar " << whichSidecar << " received the quit signal." << std::endl;
	    }
            finished = true;
	  }
	  break;

	  default:
	  {
            if(ParallelDescriptor::IOProcessor()) {
              std::cout << "**** Sidecar " << whichSidecar << " received bad signal = "
	                << sidecarSignal << std::endl;
	    }
	  }
	  break;
	}

    }
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "===== SIDECAR " << whichSidecar << " DONE. EXITING ... =====" << std::endl;
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
    int nProcs(ParallelDescriptor::NProcs());
    int ioProcNum(ParallelDescriptor::IOProcessorNumber());
    int nSidecars(0), sidecarSignal(S_SendBoxArray);
    int maxGrid(32), maxSize(16);
    int ts(0), nSteps(5);
    bool useSequential(false);
    ParmParse pp;

    if(nProcs < 8) {
      amrex::Abort("**** Error:  this test must be run with at least 8 processes.");
    }

    pp.query("maxGrid", maxGrid);
    pp.query("nComp", nComp);
    pp.query("nGhost", nGhost);
    pp.query("maxSize", maxSize);
    pp.query("nSteps", nSteps);
    pp.query("useSequential", useSequential);

    nSidecars = 3;

    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "nSidecars = " << nSidecars << std::endl;
    }

    if(ParallelDescriptor::IOProcessor()) {
      std::cout << myProcAll << ":: Resizing sidecars = " << nSidecars << std::endl;
    }

    Vector<int> randomRanks;
    if(ParallelDescriptor::IOProcessor()) {
      bool printSet(true);
      amrex::UniqueRandomSubset(randomRanks, nProcs, nProcs, printSet);
      for(int i(0); i < randomRanks.size(); ++i) {
        if(randomRanks[i] == 0) {  // ---- comprank[0] must be 0
	  randomRanks[i] = randomRanks[0];
	  randomRanks[0] = 0;
	}
      }
    }
    amrex::BroadcastArray(randomRanks, myProcAll, ioProcNum, ParallelDescriptor::Communicator());

    if(useSequential) {
      for(int i(0); i < randomRanks.size(); ++i) {
	  randomRanks[i] = i;
      }
    }

    int totalSidecarProcs(6);
    Vector<int> compProcsInAll;

    for(int i(0); i < nProcs - totalSidecarProcs; ++i) {
      compProcsInAll.push_back(randomRanks[i]);
    }

    int sCount(nProcs - totalSidecarProcs);
    Vector<Vector<int> > sidecarProcsInAll(nSidecars);
    sidecarProcsInAll[0].push_back(randomRanks[sCount++]);
    sidecarProcsInAll[0].push_back(randomRanks[sCount++]);
    sidecarProcsInAll[0].push_back(randomRanks[sCount++]);

    sidecarProcsInAll[1].push_back(randomRanks[sCount++]);

    sidecarProcsInAll[2].push_back(randomRanks[sCount++]);
    sidecarProcsInAll[2].push_back(randomRanks[sCount++]);



    bool printRanks(false);
    ParallelDescriptor::SetNProcsSidecars(compProcsInAll, sidecarProcsInAll, printRanks);

    MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;


    if(ParallelDescriptor::InSidecarGroup()) {

      int inWhichSidecar(ParallelDescriptor::InWhichSidecar());
      switch(inWhichSidecar) {
	case 0:
          SidecarEventLoop0();
	break;
	case 1:
          SidecarEventLoop1();
	break;
	case 2:
          SidecarEventLoop2();
	break;
	default:
	break;
      }

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

	if(nSidecars > 0) {

	  int whichSidecar(2);

	  if((i - ts) == 0) {  // ---- do a simple mf copy test
	    MultiFab mf(ba, comp_DM, nComp, nGhost);
	    for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
	      for(int i(0); i < mf[mfi].nComp(); ++i) {
	        mf[mfi].setVal(myProcAll + (Real) i / 1000.0, i);
	      }
	    }

	    sidecarSignal = S_CopyFabArray;
	    ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
	                              ParallelDescriptor::CommunicatorInter(whichSidecar));
	    BoxArray::SendBoxArray(ba, whichSidecar);

	    MultiFab *mfSource = &mf;
	    MultiFab *mfDest = 0;
	    bool fromComp(true);
            CopyFabArray(mfSource, mfDest, 0, 0, nComp, 0, 0, fromComp, whichSidecar);
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
	                              ParallelDescriptor::CommunicatorInter(whichSidecar));
	    BoxArray::SendBoxArray(ba, whichSidecar);  // ---- send the sidecar the unshrunk boxarray

	    MultiFab *mfSource = &mf;
	    MultiFab *mfDest = 0;
	    bool fromComp(true);
            CopyFabArray(mfSource, mfDest, 0, 0, nComp, 0, 0, fromComp, whichSidecar);
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
	                              ParallelDescriptor::CommunicatorInter(whichSidecar));
	    BoxArray::SendBoxArray(ba, whichSidecar);  // ---- send the sidecar the unshifted boxarray

	    MultiFab *mfSource = &mf;
	    MultiFab *mfDest = 0;
	    bool fromComp(true);
            CopyFabArray(mfSource, mfDest, 0, 0, nComp, 0, 0, fromComp, whichSidecar);
	  }


	  if((i - ts) == 3) {  // ---- copy part of a FabArray from the sidecar
	    MultiFab mf(ba, comp_DM, nComp, nGhost);
	    mf.setVal(-2.0);

	    sidecarSignal = S_CopyFabArrayFromSidecar;
	    ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
	                              ParallelDescriptor::CommunicatorInter(whichSidecar));

	    MultiFab *mfSource = 0;
	    MultiFab *mfDest = &mf;
	    int srcComp(4), destComp(2), numComp(1);
	    bool fromComp(false);
            CopyFabArray(mfSource, mfDest, srcComp, destComp, numComp, 0, 0, fromComp, whichSidecar);
	  }
	}
      }
      ts += nSteps;

      if(nSidecars > 0) {
	// ---- stop the sidecars
	sidecarSignal = S_QuitSignal;
	ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
	                          ParallelDescriptor::CommunicatorInter(0));
	ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
	                          ParallelDescriptor::CommunicatorInter(1));
	ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank,
	                          ParallelDescriptor::CommunicatorInter(2));
      }
    }
    
    ParallelDescriptor::Barrier();
    amrex::USleep(myProcAll/10.0);
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << myProcAll << ":: Finished timesteps" << std::endl;
    }

    ParallelDescriptor::Barrier();
    nSidecars = 0;
    ParallelDescriptor::SetNProcsSidecars(nSidecars);
    ParallelDescriptor::Barrier();


    amrex::Finalize();
    return 0;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
