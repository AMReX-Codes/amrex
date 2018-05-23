#include <iostream>
#include <AMReX.H>
#include <AMReX_Print.H>

#include <AMReX_ArrayLim.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_BaseFab.H>
#include <AMReX_BaseFab_f.H>

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
void kernel_BaseFab(amrex::BaseFab<amrex::Real> *bf1, int *val)
{
   amrex::Real value = 20.0;
   int ncomp = 1;
   val[0] = 14;
// bf1->getVal(&val[0],amrex::IntVect(0,0,0) ,0,1);
//   amrex::Box bx(amrex::IntVect(0,0,0), amrex::IntVect(13,13,13));
//   bf1->performSetVal(value, bx, 0, 1);
 
//   amrex_fort_fab_setval(bx.loVect(), bx.hiVect(), bf1->dataPtr(0), bx.loVect(), bx.hiVect(), &ncomp, &value);

//   bf1->setVal(20, bx, 0, 1);
//   bf1->getVal(&val[0],amrex::IntVect(0,0,0),0,1);

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
    // IntVect

    {
      amrex::IntVect host(1,1,1);
      amrex::IntVect *device1 = new amrex::IntVect(1,1,1);
      amrex::IntVect *device2 = new amrex::IntVect(2,2,2);

      int blockSize = 1;
      int numBlocks = 1;
//      SIMPLE_LAUNCH(numBlocks, blockSize, kernel_IntVect, device1, device2);
      kernel_IntVect<<<numBlocks,blockSize>>>(device1, device2);
      cudaDeviceSynchronize();

      amrex::Print() << "IntVect Test" << std::endl <<
                        "host = "      << host << std::endl << 
                        "device1 = "   << *device1 << std::endl <<
                        "device2 = "   << *device2 << std::endl << std::endl;
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
//      int val_gpu = 101;
//      int *d_val_gpu;
//      cudaMalloc(&d_val_gpu, size_t(sizeof(int)));
//      cudaMemcpy(d_val_gpu, &val_gpu, sizeof(int), cudaMemcpyHostToDevice);

      amrex::Real val_old = 10, val_new = 11;
      amrex::Box baseBox(amrex::IntVect(0,0,0), amrex::IntVect(13,13,13));
//      amrex::BaseFab<amrex::Real> *d_bf1 = new amrex::BaseFab<amrex::Real>(baseBox);
      amrex::BaseFab<amrex::Real> d_bf1(baseBox);

      amrex::BaseFab<int> *d_bf2 = new amrex::BaseFab<int>(baseBox);
      d_bf1.setVal(11, baseBox, 0, 1);
      d_bf1.getVal(&val_old,amrex::IntVect(1,1,1));

      amrex::Print() << "Box = " << baseBox << std::endl;

//      amrex::Real value_new = 20.0;
//      amrex_fort_fab_setval(baseBox.loVect(), baseBox.hiVect(), d_bf1->dataPtr(0), baseBox.loVect(), baseBox.hiVect(), 0, &value_new);

      int blockSize = 1;
      int numBlocks = 1;
//      kernel_BaseFab<<<numBlocks,blockSize>>>(d_bf1, d_val_gpu);
//      cudaDeviceSynchronize();
//      cudaMemcpy(&val_gpu, d_val_gpu, sizeof(int), cudaMemcpyDeviceToHost); 

 //     d_bf1->getVal(&val_new,amrex::IntVect(0,0,0));
      amrex::Print() << "BaseFab Test "  << std::endl <<
                        "pre setVal = "  << val_old << std::endl << 
//                        "gpu Val = "     << val_gpu << std::endl <<
                        "post setVal = " << val_new << std::endl; 
    }


    amrex::Finalize();
}
