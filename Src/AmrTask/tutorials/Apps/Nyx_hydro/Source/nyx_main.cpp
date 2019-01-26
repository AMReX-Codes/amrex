
#include <iostream>
#include <iomanip>
#include <sstream>

#ifndef WIN32
#include <unistd.h>
#endif

#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#ifdef BL_USE_MPI
#include <MemInfo.H>
#endif
#include <Nyx.H>
#ifdef USE_PERILLA
#include <PerillaRts.H>
#endif

#ifdef REEBER
#include <ReeberAnalysis.H> // This actually works both in situ and in-transit.
#endif

#include "Nyx_output.H"

std::string inputs_name = "";

#ifdef GIMLET
#include <DoGimletAnalysis.H>
#include <postprocess_tau_fields.H>
#include <fftw3-mpi.h>
#include <MakeFFTWBoxes.H>
#endif

#ifdef HENSON
#include <henson/context.h>
#include <henson/data.h>
#endif

using namespace amrex;

const int NyxHaloFinderSignal(42);
const int resizeSignal(43);
const int GimletSignal(55);
const int quitSignal(-44);

void
nyx_main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    // save the inputs file name for later
    if (argc > 1) {
      if (!strchr(argv[1], '=')) {
        inputs_name = argv[1];
      }
    }
    BL_PROFILE_REGION_START("main()");
    BL_PROFILE_VAR("main()", pmain);

    //
    // Don't start timing until all CPUs are ready to go.
    //
    ParallelDescriptor::Barrier("Starting main.");

    BL_COMM_PROFILE_NAMETAG("main TOP");

    Real dRunTime1 = ParallelDescriptor::second();

    std::cout << std::setprecision(10);

    int max_step;
    Real strt_time;
    Real stop_time;
    ParmParse pp;

    max_step  = -1;
    strt_time =  0.0;
    stop_time = -1.0;

    pp.query("max_step",  max_step);
    pp.query("strt_time", strt_time);
    pp.query("stop_time", stop_time);

    int how(-1);
    pp.query("how",how);

    if (strt_time < 0.0)
    {
        amrex::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0.0)
    {
        amrex::Abort("**** Error: either max_step or stop_time has to be positive!");
    }

    // Reeber has to do some initialization.
#ifdef REEBER
    reeber_int = initReeberAnalysis();
#endif

    Amr *amrptr = new Amr;
    amrptr->init(strt_time,stop_time);

#if BL_USE_MPI
    // ---- initialize nyx memory monitoring
    MemInfo *mInfo = MemInfo::GetInstance();
    mInfo->LogSummary("MemInit  ");
#endif

    const Real time_before_main_loop = ParallelDescriptor::second();

#ifdef AMREX_USE_CVODE
    Nyx::alloc_simd_vec();
#endif

    RTS rts;
    rts.Init(ParallelDescriptor::MyProc(), ParallelDescriptor::NProcs());
    rts.Iterate(amrptr, max_step, stop_time); //run coarseTimeSteps on the amr object 
    rts.Barrier();

#if 0
    bool finished(false);

    while ( ! finished) 
    {
     // If we set the regrid_on_restart flag and if we are *not* going to take
     // a time step then we want to go ahead and regrid here.
     //
     if (amrptr->RegridOnRestart()) {
       if (    (amrptr->levelSteps(0) >= max_step ) ||
               ( (stop_time >= 0.0) &&
                 (amrptr->cumTime() >= stop_time)  )    )
       {
           // Regrid only!
           amrptr->RegridOnly(amrptr->cumTime());
       }
     }

     if (amrptr->okToContinue()
          && (amrptr->levelSteps(0) < max_step || max_step < 0)
          && (amrptr->cumTime() < stop_time || stop_time < 0.0))

     {
       amrptr->coarseTimeStep(stop_time);          // ---- Do a timestep.
#ifdef HENSON
       henson_save_pointer("amr",  amrptr);        // redundant to do every timesetp, but negligible overhead
       henson_save_pointer("dmpc", Nyx::theDMPC());
       henson_yield();
#endif
     } else {
       finished = true;
     }

    }  // ---- end while( ! finished)
#endif

#ifdef AMREX_USE_CVODE
    Nyx::dealloc_simd_vec();
#endif

    const Real time_without_init = ParallelDescriptor::second() - time_before_main_loop;
    if (ParallelDescriptor::IOProcessor()) std::cout << "Time w/o init: " << time_without_init << std::endl;

    // Write final checkpoint and plotfile
    if (amrptr->stepOfLastCheckPoint() < amrptr->levelSteps(0)) {
        amrptr->checkPoint();
    }
    if (amrptr->stepOfLastPlotFile() < amrptr->levelSteps(0)) {
        amrptr->writePlotFile();
    }

    delete amrptr;

    //
    // This MUST follow the above delete as ~Amr() may dump files to disk.
    //
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    Real dRunTime2 = ParallelDescriptor::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2, IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Run time = " << dRunTime2 << std::endl;
    }

    BL_PROFILE_VAR_STOP(pmain);
    BL_PROFILE_REGION_STOP("main()");
    BL_PROFILE_SET_RUN_TIME(dRunTime2);

    amrex::Finalize();
}
