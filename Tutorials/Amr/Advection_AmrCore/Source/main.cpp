
#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include <AmrAdv.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

    // wallclock time
    const Real strt_total = ParallelDescriptor::second();

    {
	// declare an AmrAdv object to manage multilevel data
        AmrAdv amradv;
	
        // initialize AMR data
	amradv.InitData();

        // advance solution to final time
	amradv.Evolve();
	
        // wallclock time
	Real end_total = ParallelDescriptor::second() - strt_total;
	
        // print wallclock time
	ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
	if (amradv.Verbose() && ParallelDescriptor::IOProcessor()) {
	    std::cout << "\nTotal Time                     : " << end_total << '\n';
	}
    }

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
}
