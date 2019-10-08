#include <new>
#include <iostream>
#include <iomanip>

#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>
#include <PerillaRts.H>

using namespace amrex;
using namespace perilla;
int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);
    Real dRunTime1 = ParallelDescriptor::second();

    std::cout << std::setprecision(10);

    int  max_step;
    Real strt_time;
    Real stop_time;
    ParmParse pp; 

    max_step  = -1;
    strt_time =  0.0;
    stop_time = -1.0;

    pp.query("max_step",max_step);
    pp.query("strt_time",strt_time);
    pp.query("stop_time",stop_time);

    if (strt_time < 0.0) {
        amrex::Abort("MUST SPECIFY a non-negative strt_time"); 
    }

    if (max_step < 0 && stop_time < 0.0) {
	amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    {
        Amr amr;

        amr.init(strt_time,stop_time);

#if USE_PERILLA_PTHREADS
        RTS rts;
        rts.Init(ParallelDescriptor::MyProc(), ParallelDescriptor::NProcs());
        rts.Iterate(&amr, max_step, stop_time); //run coarseTimeSteps on the amr object 

        rts.Barrier();
        // Write final checkpoint and plotfile
        if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
            amr.checkPoint();
        }

        if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
            amr.writePlotFile();
        }
        rts.Finalize();
#elif defined(USE_PERILLA_ON_DEMAND)
        while ( amr.okToContinue() &&
               (amr.levelSteps(0) < max_step || max_step < 0) &&
               (amr.cumTime() < stop_time || stop_time < 0.0) )

        {
            //
            // Do a coarse timestep.  Recursively calls timeStep()
            //
            amr.coarseTimeStep(stop_time);
        }

        // Write final checkpoint and plotfile
        if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
            amr.checkPoint();
        }

        if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
            amr.writePlotFile();
        }
#endif
    }

    Real dRunTime2 = ParallelDescriptor::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Run time = " << dRunTime2 << std::endl;
    }

    amrex::Finalize();

    return 0;
}
