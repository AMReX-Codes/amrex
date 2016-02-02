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

#include <Geometry.H>
#include <ParallelDescriptor.H>
#include <ParmParse.H>
#include <RealBox.H>
#include <Utility.H>
#include <ParmParse.H>

#include <BoxLib.H>
#include <InTransitAnalysis.H>

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
    BL_ASSERT(MySignal != ParallelDescriptor::SidecarQuitSignal);
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
      Geometry::SendGeometryToSidecars(&geom);
      ParallelDescriptor::Bcast(&time_step, 1, 0, ParallelDescriptor::CommunicatorInter());

      InTransitAnalysis ita;

      ita.Initialize(mf, geom, time_step);
      ita.DoAnalysis();
      ita.Finalize();

      if (ParallelDescriptor::IOProcessor())
        std::cout << "Sidecars completed analysis." << std::endl;
    }
    else if (in_signal != ParallelDescriptor::SidecarQuitSignal)
    {
      std::ostringstream ss_error_msg;
      ss_error_msg << "Unknown signal sent to sidecars: -----> " << signal << " <-----" << std::endl;
      BoxLib::Error(const_cast<const char*>(ss_error_msg.str().c_str()));
    }

    return out_signal;
  }

  static void STATIC_INIT () {
    if (ParallelDescriptor::InSidecarGroup() && ParallelDescriptor::IOProcessor())
        std::cout << "This is where the signal handler would initialize a bunch of stuff ..." << std::endl;
    ParallelDescriptor::AddSignalHandler(STATIC_SIGNAL_HANDLER);
  }

  static void STATIC_CLEAN () {
    if (ParallelDescriptor::InSidecarGroup() && ParallelDescriptor::IOProcessor())
        std::cout << "This is where the signal handler would clean stuff up ..." << std::endl;
  }

  static void RunAtStatic () {
    BoxLib::ExecOnInitialize(STATIC_INIT);
    BoxLib::ExecOnFinalize(STATIC_CLEAN);
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
    const int nSidecarProcs(2);
    ParallelDescriptor::SetNProcsSidecar(nSidecarProcs);
#endif

    BoxLib::Initialize(argc,argv);

    // The sidecar group has already called BoxLib::Finalize() by the time we
    // are out of BoxLib::Initialize(), so make them quit immediately.
    // Everything below this point is done on the compute group only.
#ifdef IN_TRANSIT
    if (ParallelDescriptor::InSidecarGroup())
        return 0;
#endif

    // A flag you need for broadcasting across MPI groups. We always broadcast
    // the data to the sidecar group from the IOProcessor on the compute group.
#ifdef IN_TRANSIT
    const int MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;
#endif

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
    MultiFab mf(ba, nComp, nGhost, comp_DM, Fab_allocate);
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
    int signal;

    // Pretend we're doing a halo-finding analysis for Nyx.
    signal = MySignal;

    // Pretend we're looping over time steps.
    for (unsigned int i = 0; i < 20; ++i)
    {
        ParallelDescriptor::Bcast(&signal, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
        // For this particular analysis we need to send a MultiFab, a Geometry,
        // and a time step index to the sidecars. For other analyses you will
        // need to send different data.
        MultiFab::SendMultiFabToSidecars(&mf);
        Geometry::SendGeometryToSidecars(&geom);
        ParallelDescriptor::Bcast(&i, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
    }

    // Don't forget to tell the sidecars to quit! Otherwise they'll keep
    // waiting for signals for more work to do.
    signal = ParallelDescriptor::SidecarQuitSignal;
    ParallelDescriptor::Bcast(&signal, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
#endif

    BoxLib::Finalize();
    return 0;
}
