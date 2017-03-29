
#include <iostream>
#include <mpi.h>

#include <WarpXWrappers.h>

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    std::cout << "MPI is initialized " << std::endl;

    amrex_init_with_inited_mpi(argc, argv, MPI_COMM_WORLD);

    std::cout << "AMReX is initialized" << std::endl;

    warpx_init();

    warpx_evolve(-1);

    warpx_finalize();

    amrex_finalize(0);

    MPI_Finalize();
}
