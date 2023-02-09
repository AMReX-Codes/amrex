#include <AMReX.H>


int main (int argc, char* argv[])
{
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif
    {
        amrex::Initialize(argc,argv);
        amrex::Finalize();
    }

    {
        amrex::Initialize(argc,argv);
        amrex::Finalize();
    }

    {
        amrex::Initialize(argc,argv);
        amrex::Finalize();
    }
#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif
}
