#include <iostream>
#include <AMReX.H>
#include <AMReX_Print.H>

#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>

DEVICE_FUNCTION
void kernel_Vector()
{
}


DEVICE_FUNCTION
void kernel_IntVect(amrex::IntVect *iv1, amrex::IntVect *iv2)
{
   *iv2 = *iv1 + *iv2;

   amrex::IntVect local(3, 3, 3);
   *iv1 = *iv1 + local; 
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    // Simple cuda action to make sure all tests have cuda.
    // Allows nvprof to return data.
    int devices;
    cudaGetDeviceCount(&devices);
    amrex::Print() << "Hello world from AMReX version " << amrex::Version() << ". GPU devices: " << devices << "\n";
    // ===================================
    /*
      Priority List
      -------------
      0) AMReX_Vector
      1) IntVect
      2) Box 
      3) BaseFab
      4) FabArray
      5) Geometry 
    */
    // ===================================
    // AMReX_Vector
    {
      amrex::Vector *vector1 = new amrex::vector(10, -1);

      int blockSize = 256;
      int numBlocks = 1;
      kernel_Vector<<<numBlocks,blockSize>>>(device1, device2);
      cudaDeviceSynchronize();
    }

    // ===================================
    // IntVect
    {
      amrex::IntVect host(1,1,1);
      amrex::IntVect *device1 = new amrex::IntVect(1,1,1);
      amrex::IntVect *device2 = new amrex::IntVect(2,2,2);

      int blockSize = 256;
      int numBlocks = 1;
//      SIMPLE_LAUNCH(numBlocks, blockSize, kernel_IntVect, device1, device2);
      kernel_IntVect<<<numBlocks,blockSize>>>(device1, device2);
      cudaDeviceSynchronize();

      amrex::Print() << "host = " << host << std::endl << 
                        "device1 = " << *device1 << std::endl <<
                        "device2 = " << *device2 << std::endl;
    }



    amrex::Finalize();
}
