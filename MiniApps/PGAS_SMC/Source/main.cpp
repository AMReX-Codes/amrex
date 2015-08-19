#include <iostream>
#include <iomanip>

#include <unistd.h>

#include <REAL.H>
#include <Utility.H>
#include <IntVect.H>
#include <Box.H>
#include <Amr.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <AmrLevel.H>

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    Real dRunTime1 = ParallelDescriptor::second();

    std::cout << std::setprecision(10);

    int  max_step       = -1;
    Real strt_time      = 0.0;
    Real stop_time      = -1.0;

    ParmParse pp; 

    pp.get("max_step",  max_step);
    pp.get("strt_time", strt_time);
    pp.get("stop_time", stop_time);

    if (strt_time < 0.0)
    {
        BoxLib::Abort("MUST SPECIFY a non-negative strt_time"); 
    }

    if (max_step < 0 && stop_time < 0.0) {
	BoxLib::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    Amr* amrptr = new Amr;

    amrptr->init(strt_time,stop_time);

    Real dRunTime2 = ParallelDescriptor::second();

    while ( amrptr->okToContinue()                            &&
	    (amrptr->levelSteps(0) < max_step || max_step < 0) &&
	    (amrptr->cumTime() < stop_time || stop_time < 0.0) )
    {
        //
        // Do a timestep.
        //
        amrptr->coarseTimeStep(stop_time);
    }

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

    Real dRunTime3 = ParallelDescriptor::second();

    Real dRT[2] = {dRunTime3-dRunTime1, dRunTime3-dRunTime2};

    ParallelDescriptor::ReduceRealMax(dRT,2,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Run time          = " << dRT[0] << std::endl;
        std::cout << "Run time w/o init = " << dRT[1] << std::endl;
    }

    if (CArena* arena = dynamic_cast<CArena*>(BoxLib::The_Arena()))
    {
        //
        // A barrier to make sure our output follows that of RunStats.
        //
        ParallelDescriptor::Barrier();
        //
        // We're using a CArena -- output some FAB memory stats.
        //
        // This'll output total # of bytes of heap space in the Arena.
        //
        // It's actually the high water mark of heap space required by FABs.
        //
        char buf[256];

        sprintf(buf,
                "CPU(%d): Heap Space (bytes) used by Coalescing FAB Arena: %ld",
                ParallelDescriptor::MyProc(),
                arena->heap_space_used());

        std::cout << buf << std::endl;
    }

    BL_PROFILE_VAR_STOP(pmain);

    BoxLib::Finalize();

    return 0;
}
