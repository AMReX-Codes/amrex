
#include <new>
#include <iostream>
#include <iomanip>

#include <Amr.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <AmrLevel.H>

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);
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
        BoxLib::Abort("MUST SPECIFY a non-negative strt_time"); 
    }

    if (max_step < 0 && stop_time < 0.0) {
	BoxLib::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    Amr* amrptr = new Amr;

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

    // Write final checkpoint and plotfile
    if (amrptr->stepOfLastCheckPoint() < amrptr->levelSteps(0)) {
        amrptr->checkPoint();
    }

    if (amrptr->stepOfLastPlotFile() < amrptr->levelSteps(0)) {
        amrptr->writePlotFile();
    }

    delete amrptr;

    Real dRunTime2 = ParallelDescriptor::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Run time = " << dRunTime2 << std::endl;
    }

    BoxLib::Finalize();

    return 0;
}
