
#include <iostream>

#include <BoxLib.H>
#include <BLProfiler.H>
#include <ParallelDescriptor.H>

#include <AmrAdv.H>

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    const Real strt_total = ParallelDescriptor::second();

    {
	AmrAdv amradv;
	
	amradv.InitData();

	amradv.Evolve();
	
	Real end_total = ParallelDescriptor::second() - strt_total;
	
	ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
	if (amradv.Verbose() && ParallelDescriptor::IOProcessor()) {
	    std::cout << "\nTotal Time                     : " << end_total << '\n';
	}
    }

    BL_PROFILE_VAR_STOP(pmain);

    BoxLib::Finalize();
}
