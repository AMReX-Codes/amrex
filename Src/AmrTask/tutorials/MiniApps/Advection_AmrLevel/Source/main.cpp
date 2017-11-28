
#include <new>
#include <iostream>
#include <iomanip>

#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>
#include "AMReX_AsyncMFIter.H"
#include <AMReX_AmrTask.H>

using namespace amrex;

class MyAction :public AmrTask{
    public:
	virtual Real advanceTask (Real time,
		Real dt,
		int  iteration,
		int  ncycle){}
};


int main (int   argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    Real dRunTime1 = ParallelDescriptor::second();

    int  max_step;
    Real strt_time;
    Real stop_time;

    {
	ParmParse pp; 

	max_step  = -1;
	strt_time =  0.0;
	stop_time = -1.0;

	pp.query("max_step",max_step);
	pp.query("strt_time",strt_time);
	pp.query("stop_time",stop_time);
    }

    if (strt_time < 0.0) {
	amrex::Abort("MUST SPECIFY a non-negative strt_time"); 
    }

    if (max_step < 0 && stop_time < 0.0) {
	amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    {
	Amr amr;

	amr.init(strt_time,stop_time);

	AMFIter<MyAction> mfi(&amr, max_step, stop_time, false/*do tiling*/);
	mfi.Iterate();

#if 0
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

    amrex::Print() << "Run time = " << dRunTime2 << std::endl;

    amrex::Finalize();

    return 0;
}
