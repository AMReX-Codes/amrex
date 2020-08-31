
#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include <AmrCoreAdv.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
        // timer for profiling
        BL_PROFILE("main()");

        // wallclock time
        const Real strt_total = amrex::second();

        // constructor - reads in parameters from inputs file
        //             - sizes multilevel arrays and data structures
        AmrCoreAdv amr_core_adv;
	
        // initialize AMR data
	amr_core_adv.InitData();

        // advance solution to final time
	amr_core_adv.Evolve();
	
        // wallclock time
	Real end_total = amrex::second() - strt_total;
	
	if (amr_core_adv.Verbose()) {
            // print wallclock time
            ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
            amrex::Print() << "\nTotal Time: " << end_total << '\n';
	}
    }

    amrex::Finalize();
}
