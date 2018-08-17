#include <AMReX.H>
#include <AMReX_Geometry.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    
     Geometry geom;
    
    amrex::Finalize();
}
