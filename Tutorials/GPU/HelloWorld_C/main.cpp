#include <iostream>
#include <AMReX.H>
#include <AMReX_Print.H>

#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_BaseFab.H>

DEVICE_FUNCTION
void kernel_Vector(amrex::Vector<int> vec1, amrex::Vector<int> vec2)
{
/*   for (int i=0; i<vec1.size(); ++i)
   {
     vec1[i] = vec1[i] + vec2[i];
   }
*/
}

DEVICE_FUNCTION
void kernel_IntVect(amrex::IntVect *iv1, amrex::IntVect *iv2)
{
   *iv2 = *iv1 + *iv2;

   amrex::IntVect local(3, 3, 3);
   *iv1 = *iv1 + local; 
}

DEVICE_FUNCTION
void kernel_Box(amrex::Box *box1, amrex::Box *box2)
{

}

DEVICE_FUNCTION
void kernel_BaseFab(amrex::BaseFab<int> *bf1, amrex::BaseFab<int> *bf2)
{
   amrex::Box baseBox(amrex::IntVect(0,0,0), amrex::IntVect(9,9,9));
   bf1->setVal(1, baseBox, 1, 1);
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
      0-1) Functions with std calls.
      1) IntVect
      2) IndexType
      3) Box 
      4) BaseFab
      5) FabArray
      6) Geometry 

      Top Down
      --------
      A) Setval
      B) Copy
      C) Fortran Functions
    */
    // ===================================
    // AMReX_Vector
    //{
    //  amrex::Vector<int> *vector1 = new amrex::Vector<int>(10, -1);
    //  amrex::Vector<int> *vector2 = new amrex::Vector<int>(10, 5);
    //
    //  int blockSize = 256;
    //  int numBlocks = 1;
    //  kernel_Vector<<<numBlocks,blockSize>>>(*vector1, *vector2);
    //  cudaDeviceSynchronize();
    //
    //  amrex::Print() << "vector1 :: vector 2:" << std::endl;
    //  for (int i=0; i<vector1->size(); ++i)
    //  {
    //    amrex::Print() << (*vector1)[i] << " :: " << (*vector2)[i] << std::endl;
    //  }
    //}
    // ===================================
    // Array 
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
/*
    // ===================================
    // Box 
    {
      amrex::Box host_box(amrex::IntVect(0,0,0), amrex::IntVect(9,9,9));
      amrex::Box *d_box1 = new amrex::Box(amrex::IntVect(1,1,1), amrex::IntVect(5,5,5));
      amrex::Box *d_box2 = new amrex::Box(amrex::IntVect(3,3,3), amrex::IntVect(6,6,6));

      int blockSize = 256;
      int numBlocks = 1;
      kernel_Box<<<numBlocks,blockSize>>>(d_box1, d_box2);
      cudaDeviceSynchronize();

      amrex::Print() << "host = " << host << std::endl << 
                        "device1 = " << *d_box1 << std::endl <<
                        "device2 = " << *d_box2 << std::endl;
    }
*/
    // ===================================
    // BaseFab
    {
      amrex::Box baseBox(amrex::IntVect(0,0,0), amrex::IntVect(13,13,13));
      amrex::BaseFab<int> *d_bf1 = new amrex::BaseFab<int>(baseBox);
      amrex::BaseFab<int> *d_bf2 = new amrex::BaseFab<int>(baseBox);
//      amrex::BaseFab<int> h_bf(baseBox);

      int blockSize = 256;
      int numBlocks = 1;
      kernel_BaseFab<<<numBlocks,blockSize>>>(d_bf1, d_bf2);
      cudaDeviceSynchronize();
    }


    amrex::Finalize();
}
