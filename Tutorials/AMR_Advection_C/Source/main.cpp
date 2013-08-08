//
// $Id: main.cpp,v 1.1 2010-11-18 22:29:21 almgren Exp $
//
#include <winstd.H>

#include <new>
#include <cstdio>
#include <iostream>
#include <iomanip>

#ifndef WIN32
#include <unistd.h>
#endif

#include <CArena.H>
#include <REAL.H>
#include <Utility.H>
#include <IntVect.H>
#include <Box.H>
#include <Amr.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <AmrLevel.H>

#ifdef USE_SDCLIB
#include <MLSDCAmr.H>
#endif

int
main (int   argc,
      char* argv[])
{

    //
    // Make sure to catch new failures.
    //
    BoxLib::Initialize(argc,argv);
    Real dRunTime1 = ParallelDescriptor::second();

#ifdef BL_XT3
    const int siobsize(8192);
    char stdiobufff[siobsize];
    std::cout.rdbuf()->pubsetbuf(stdiobufff, siobsize);
#endif

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

    if (strt_time < 0.0)
    {
        BoxLib::Abort("MUST SPECIFY a non-negative strt_time"); 
    }

    if (max_step < 0 && stop_time < 0.0) {
      BoxLib::Abort(
       "Exiting because neither max_step nor stop_time is non-negative.");
    }

#ifdef USE_SDCLIB
    Amr* amrptr = new MLSDCAmr;
#else
    Amr* amrptr = new Amr;
#endif

    amrptr->init(strt_time,stop_time);

    // If we set the regrid_on_restart flag and if we are *not* going to take
    //    a time step then we want to go ahead and regrid here.
    if ( amrptr->RegridOnRestart() && 
         ( (amrptr->levelSteps(0) >= max_step) ||
           (amrptr->cumTime() >= stop_time) ) )
           {
           //
           // Regrid only!
           //
           amrptr->RegridOnly(amrptr->cumTime());
           }

    while ( amrptr->okToContinue()                            &&
           (amrptr->levelSteps(0) < max_step || max_step < 0) &&
           (amrptr->cumTime() < stop_time || stop_time < 0.0) )

    {
        //
        // Do a timestep.
        //
        amrptr->coarseTimeStep(stop_time);
    }

    delete amrptr;
    //
    // This MUST follow the above delete as ~Amr() may dump files to disk.
    //
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    Real dRunTime2 = ParallelDescriptor::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Run time = " << dRunTime2 << std::endl;
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

    BoxLib::Finalize();

    return 0;
}
