
#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include <AmrAdv.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

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

    amrex::Finalize();
}
