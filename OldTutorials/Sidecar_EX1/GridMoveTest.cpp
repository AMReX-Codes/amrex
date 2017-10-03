// An example showing how to clone a MultiFab from the compute MPI group to the
// "sidecar" group and do analysis on the MultiFab on the fly.

#include <new>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sstream>
#include <unistd.h>

#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_RealBox.H>
#include <AMReX_Utility.H>
#include <AMReX_ParmParse.H>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <InTransitAnalysis.H>

using namespace amrex;

const int SidecarQuitSignal = -51;

// In this anonymous namespace we define the workflow which occurs when the
// sidecars receive a particular (user-defined) signal. You have a lot of
// flexibility here. You could define all of your workflow directly in this
// namespace, which may be useful if you have only simple operations to
// perform. Or you could define it all within functions (or member functions of
// a class), so that this namespace would consist of function calls and little
// else. This tutorial presents a mixture of both.

namespace
{
#ifdef IN_TRANSIT
  const int MySignal = 42;

  static int STATIC_SIGNAL_HANDLER (int in_signal) {
    BL_ASSERT(MySignal != SidecarQuitSignal);
    int out_signal = in_signal;
    if (in_signal == MySignal)
    {
      if (ParallelDescriptor::IOProcessor())
        std::cout << "Sidecars got the analysis signal!" << std::endl;

      MultiFab mf;
      Geometry geom;
      int time_step;
      // Receive the necessary data for doing analysis.
      MultiFab::SendMultiFabToSidecars(&mf);
      Geometry::SendGeometryToSidecar(&geom, 0);
      ParallelDescriptor::Bcast(&time_step, 1, 0, ParallelDescriptor::CommunicatorInter(0));

/*
      MultiFab::copyInter(mfSource, mfDest, srcComp, destComp, numComp,
                        srcNGhost, destNGhost,
                        commSrc, commDest, commInter, commBoth,
                        isSrc);
*/

      InTransitAnalysis ita;

      ita.Initialize(mf, geom, time_step);
      ita.DoAnalysis();
      ita.Finalize();

      if (ParallelDescriptor::IOProcessor())
        std::cout << "Sidecars completed analysis." << std::endl;
    }
    else if (in_signal != SidecarQuitSignal)
    {
      std::ostringstream ss_error_msg;
      ss_error_msg << "Unknown signal sent to sidecars: -----> " << in_signal << " <-----" << std::endl;
      amrex::Error(const_cast<const char*>(ss_error_msg.str().c_str()));
    }

    return out_signal;
  }

  static void STATIC_INIT () {
    std::cout << "_in STATIC_INIT" << std::endl;
    if (ParallelDescriptor::InSidecarGroup() && ParallelDescriptor::IOProcessor())
        std::cout << "This is where the signal handler would initialize a bunch of stuff ..." << std::endl;
    ParallelDescriptor::AddSignalHandler(STATIC_SIGNAL_HANDLER);
  }

  static void STATIC_CLEAN () {
    std::cout << "_in STATIC_CLEAN" << std::endl;
    if (ParallelDescriptor::InSidecarGroup() && ParallelDescriptor::IOProcessor())
        std::cout << "This is where the signal handler would clean stuff up ..." << std::endl;
  }

  static void RunAtStatic () {
    amrex::ExecOnInitialize(STATIC_INIT);
    amrex::ExecOnFinalize(STATIC_CLEAN);
  };
#endif
}



// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {

  RunAtStatic();

    // Unfortunately the # of sidecars is currently a compile-time constant
    // because ParmParse must come AFTER Initialize(), but setting the # of
    // sidecars must come BEFORE.

    // TODO: change Initialize() so that we can read # of sidecars from an
    // inputs file.

#ifdef IN_TRANSIT
    //const int nSidecarProcs(2);
    //ParallelDescriptor::SetNProcsSidecar(nSidecarProcs);
#endif

    amrex::Initialize(argc,argv);

    int myProcAll(ParallelDescriptor::MyProcAll());

std::cout << myProcAll << "::_here 0." << std::endl;
#ifdef IN_TRANSIT
    const int nSidecarProcs(2);
    ParallelDescriptor::SetNProcsSidecars(nSidecarProcs);
#endif

std::cout << myProcAll << "::_here 1." << std::endl;
    // The sidecar group has already called amrex::Finalize() by the time we
    // are out of amrex::Initialize(), so make them quit immediately.
    // Everything below this point is done on the compute group only.
#ifdef IN_TRANSIT
    if (ParallelDescriptor::InSidecarGroup()) {
	std::cout << myProcAll << "::**** returning from sidecar group." << std::endl;
	sleep(4);
        return 0;
    }
#endif

std::cout << myProcAll << "::_here 2." << std::endl;
    // A flag you need for broadcasting across MPI groups. We always broadcast
    // the data to the sidecar group from the IOProcessor on the compute group.
#ifdef IN_TRANSIT
    const int MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;
#endif

std::cout << myProcAll << "::_here 3." << std::endl;
ParallelDescriptor::Barrier();
    int maxGrid;
    int nComp;
    int nGhost;
    int maxSize;

    ParmParse pp;
    pp.get("maxGrid", maxGrid);
    pp.get("nComp", nComp);
    pp.get("nGhost", nGhost);
    pp.get("maxSize", maxSize);

    // Build a Box to span the problem domain, and split it into a BoxArray
    Box baseBox(IntVect(0,0,0), IntVect(maxGrid-1, maxGrid-1, maxGrid-1));
    BoxArray ba(baseBox);
    ba.maxSize(maxSize);

    // This is the DM for the compute processes.
    DistributionMapping comp_DM;
    comp_DM.define(ba, ParallelDescriptor::NProcsComp());

    // Make a MultiFab and populate it with a bunch of random numbers.
    MultiFab mf(ba, comp_DM, nComp, nGhost);
    for(int i(0); i < mf.nComp(); ++i) {
      mf.setVal(rand()%100, i, 1);
    }

    // This defines the physical size of the box.
    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
      real_box.setLo(n, 0.0);
      real_box.setHi(n, 42.0);
    }

    // This says we are using Cartesian coordinates
    int coord(0);

    // This sets the boundary conditions to be non-periodic.
    int is_per[BL_SPACEDIM];
    for (int n = 0; n < BL_SPACEDIM; n++) is_per[n] = 0;

    // This defines a Geometry object which is useful for writing the plotfiles
    Geometry geom(baseBox, &real_box, coord, is_per);

#ifdef IN_TRANSIT
    // The signal for telling the sidecars what to do.
    int sidecarSignal;

    // Pretend we're doing a halo-finding analysis for Nyx.
    sidecarSignal = MySignal;

std::cout << myProcAll << "::_here 4." << std::endl;
ParallelDescriptor::Barrier();
    // Pretend we're looping over time steps.
    for (unsigned int i = 0; i < 20; ++i)
    {
std::cout << myProcAll << "::_here 5 i = " << i << std::endl;
ParallelDescriptor::Barrier();
        ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter(0));
        // For this particular analysis we need to send a MultiFab, a Geometry,
        // and a time step index to the sidecars. For other analyses you will
        // need to send different data.
        std::cout << myProcAll << "::_here 6 i = " << i << std::endl;
        ParallelDescriptor::Barrier();
        MultiFab::SendMultiFabToSidecars(&mf);
        std::cout << myProcAll << "::_here 7 i = " << i << std::endl;
        ParallelDescriptor::Barrier();
        Geometry::SendGeometryToSidecar(&geom, 0);
        ParallelDescriptor::Bcast(&i, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter(0));
    }

    // Don't forget to tell the sidecars to quit! Otherwise they'll keep
    // waiting for signals for more work to do.
    sidecarSignal = SidecarQuitSignal;
    ParallelDescriptor::Bcast(&sidecarSignal, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter(0));
#endif

    std::cout << "_calling Finalize()" << std::endl;
    amrex::Finalize();
    return 0;
}
