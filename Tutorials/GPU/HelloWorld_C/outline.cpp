#include <AMReX.H>
#include "MyTest.H"


__global__ void kernel_in_cpp()
{
    Box bx(IntVect(0), IntVect(15));
    BaseFab<int> fab(bx, 1);
    fab.setVal(0);
}


__global__ void kernel2_in_cpp(BaseFab<int>& fab)
{
    fab.setVal(..);
    .plus
}

for (MFIter mfi(mf); mfi.isValid(); ++mfi)
{
    kernel_in_cpp<<<...>>>(mf[mfi]);
}

{
  Box     domain;   // My index space.
  IntVect dlen;     // Length of domain in each direction.
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    // Simple cuda action for profiling.
    int devices;
    cudaGetDeviceCount(&devices);
    amrex::Print() << "Hello world from AMReX version " << amrex::Version() << ". GPU devices: " << devices << "\n";
    // ===================================
    {
        kernel_in_cpp();
    }

    {
        BaseFab<int>* fab = new BaseFab<int>(....);
        kernel2_in_cpp (*fab);

        amrex::Print() << fab << "\n";
    }

    amrex::Finalize();
}
