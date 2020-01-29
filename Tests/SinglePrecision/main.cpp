
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        static_assert(sizeof(amrex::Real) == sizeof(float), "amrex::Real and float have different size");
        amrex::Print() << "Hello world from AMReX version " << amrex::Version() << "\n";
    }
    amrex::Finalize();
}

