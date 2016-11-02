
#include <iostream>
#include <BLProfiler.H>
#include <ParallelDescriptor.H>
#include <WarpXWrappers.h>

int main(int argc, char* argv[])
{
#ifdef BOXLIB_INIT_MPI
    boxlib_init(&argc, &argv);
#else
    MPI_Init(&argc, &argv);
    boxlib_init_with_inited_mpi(&argc, &argv, MPI_COMM_WORLD);
#endif

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

#ifdef BOXLIB_INIT_MPI
    boxlib_finalize(1);
#else
    boxlib_finalize(0);
    MPI_Finalize();
#endif
}
