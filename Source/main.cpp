
#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include <WarpX.H>

using namespace amrex;

int main(int argc, char* argv[])
{
#if defined(_OPENMP) && defined(WARPX_USE_PSATD)
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    assert(provided >= MPI_THREAD_FUNNELED);
#else
    MPI_Init(&argc, &argv);
#endif

    amrex::Initialize(argc,argv,MPI_COMM_WORLD);

    BL_PROFILE_VAR("main()", pmain);

    const Real strt_total = ParallelDescriptor::second();

    {
	WarpX warpx;
	
	warpx.InitData();

	warpx.Evolve();
	
	Real end_total = ParallelDescriptor::second() - strt_total;
	
	ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
	if (warpx.Verbose()) {
            amrex::Print() << "Total Time                     : " << end_total << '\n';
	}
    }

    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
    MPI_Finalize();
}
