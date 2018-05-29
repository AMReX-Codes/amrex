#include <iostream>
#include <AMReX.H>
#include <AMReX_Print.H>

#include <AMReX_ArrayLim.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_BaseFab.H>
#include <AMReX_BaseFab_f.H>

DEVICE_FUNCTION
void kernel_IntVect(amrex::IntVect *iv1, amrex::IntVect *iv2)
{
   *iv2 = *iv1 + *iv2;

   amrex::IntVect local = *iv2;
   *iv1 = *iv1 + local; 
}

DEVICE_FUNCTION
void kernel_BaseFab(amrex::BaseFab<amrex::Real> *bf1, int *val, amrex::IntVect *iv)
{
  *iv = bf1->box().bigEnd();
  *val = bf1->nComp();
}

int main (int argc, char* argv[])
{
    std::cout << "**********************************\n";
    amrex::Initialize(argc, argv);
    amrex::Print() << "amrex::Initialize complete." << "\n";

    // ===================================
    // Simple cuda action to make sure all tests have cuda.
    // Allows nvprof to return data.
    int devices;
    cudaGetDeviceCount(&devices);
    amrex::Print() << "Hello world from AMReX version " << amrex::Version() << ". GPU devices: " << devices << "\n";
    amrex::Print() << "**********************************\n"; 
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
                        "(1,1,1) + (2,2,2) = "   << *device2 << std::endl <<
                        " + (1,1,1)  = "   << *device1 << std::endl <<
                        "**********************************" << std::endl << std::endl;
    }

    // ===================================
    // BaseFab
    {
      int val_gpu = 101;
      int *d_val_gpu;
      cudaMalloc(&d_val_gpu, size_t(sizeof(int)));
      cudaMemcpy(d_val_gpu, &val_gpu, sizeof(int), cudaMemcpyHostToDevice);

      amrex::IntVect *big = new amrex::IntVect(1,1,1);
      amrex::Real val_old = 10, val_new = 11;
      amrex::Box baseBox(amrex::IntVect(0,0,0), amrex::IntVect(14,14,14));
      amrex::BaseFab<amrex::Real> *d_bf1 = new amrex::BaseFab<amrex::Real>(baseBox);
//      amrex::BaseFab<amrex::Real> d_bf1(baseBox);

      amrex::BaseFab<int> *d_bf2 = new amrex::BaseFab<int>(baseBox);
      d_bf1->setVal(11, baseBox, 0, 1);
      d_bf1->getVal(&val_old,amrex::IntVect(1,1,1));

      amrex::Print() << "Box = " << baseBox << std::endl;

      cudaPointerAttributes ptr_attr;
      cudaPointerGetAttributes(&ptr_attr, d_bf1);
      amrex::Print() << "d_bf1.isManaged = " << ptr_attr.isManaged << std::endl;
      cudaPointerGetAttributes(&ptr_attr, d_bf1->dataPtr());
      amrex::Print() << "d_bf1.dataPtr().isManaged = " << ptr_attr.isManaged << std::endl;

//  amrex::Real value_new = 20.0;
//  amrex_fort_fab_setval(baseBox.loVect(), baseBox.hiVect(), d_bf1->dataPtr(0), baseBox.loVect(), baseBox.hiVect(), 0, &value_new);

      int blockSize = 1;
      int numBlocks = 1;
      kernel_BaseFab<<<numBlocks,blockSize>>>(d_bf1, d_val_gpu, big);
      cudaDeviceSynchronize();
      cudaMemcpy(&val_gpu, d_val_gpu, sizeof(int), cudaMemcpyDeviceToHost); 

 //     d_bf1->getVal(&val_new,amrex::IntVect(0,0,0));
      amrex::Print() << "BaseFab Test "  << std::endl <<
                        "Box::bigEnd = " << *big << std::endl <<
                        "pre setVal = "  << val_old << std::endl << 
                        "gpu Val = "     << val_gpu << std::endl <<
                        "post setVal = " << val_new << std::endl <<
                        "**********************************" << std::endl << std::endl;
    }

    amrex::Finalize();
}
