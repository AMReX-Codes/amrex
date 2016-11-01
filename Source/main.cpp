
#include <iostream>
#include <BLProfiler.H>
#include <ParallelDescriptor.H>
#include <WarpXWrappers.h>

int main(int argc, char* argv[])
{
    boxlib_init(argc, &argv);

    BL_PROFILE_VAR("main()", pmain);	

    const Real strt_total = ParallelDescriptor::second();

    warpx_init();
    warpx_evolve(-1);
    warpx_finalize();

    Real end_total = ParallelDescriptor::second() - strt_total;
    
    ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Total Time                     : " << end_total << '\n';
    }

    BL_PROFILE_VAR_STOP(pmain);

    boxlib_finalize();
}
