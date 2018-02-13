
#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include <AmrCoreAdv.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

    // wallclock time
    const Real strt_total = ParallelDescriptor::second();

    {
        // constructor - reads in parameters from inputs file
        //             - sizes multilevel arrays and data structures
        AmrCoreAdv amr_core_adv;
	
        // initialize AMR data
	amr_core_adv.InitData();

        // advance solution to final time
	amr_core_adv.Evolve();
	
        // wallclock time
	Real end_total = ParallelDescriptor::second() - strt_total;
	
        // print wallclock time
	ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
	if (amr_core_adv.Verbose()) {
            amrex::Print() << "\nTotal Time: " << end_total << '\n';
	}
    }

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
}
