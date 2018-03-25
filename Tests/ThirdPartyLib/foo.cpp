#include <mpi.h>
#include <AMReX.H>
#include <AMReX_Print.H>

extern "C"
{
    void foo(MPI_Comm comm)
    {
        amrex::Initialize(comm);
        
        amrex::Print() << " foo: AMReX C++ has been initialized.\n";
        
        amrex::Finalize();

        // After amrex::Finalize(), amrex can no longer be used.
        int rank;
        MPI_Comm_rank(comm, &rank);
        if (rank == 0) {
            std::cout << " foo: AMReX C++ has been finalized.\n";
        }
    }
}
