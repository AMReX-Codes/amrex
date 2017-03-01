
#include <iostream>

#include <BoxLib.H>
#include <BLProfiler.H>
#include <ParallelDescriptor.H>

#include <WarpX.H>

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    const Real strt_total = ParallelDescriptor::second();

    {
	WarpX warpx;
	
	warpx.InitData();

    for (int i = 0; i<100; ++i) {
	warpx.Evolve(10);
    }
	
	Real end_total = ParallelDescriptor::second() - strt_total;
	
	ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
	if (warpx.Verbose() && ParallelDescriptor::IOProcessor()) {
	    std::cout << "Total Time                     : " << end_total << '\n';
	}
    }

    BL_PROFILE_VAR_STOP(pmain);

    BoxLib::Finalize();
}
