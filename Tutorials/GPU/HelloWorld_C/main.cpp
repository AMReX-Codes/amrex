#include <iostream>
#include <AMReX.H>
#include <AMReX_Print.H>

#include <AMReX_ArrayLim.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_BaseFab.H>
#include <AMReX_BaseFab_f.H>

AMREX_CUDA_GLOBAL
void kernel_IntVect(amrex::IntVect *iv1, amrex::IntVect *iv2)
{
   *iv2 = *iv1 + *iv2;

   amrex::IntVect local = *iv2;
   *iv1 = *iv1 + local; 
}

AMREX_CUDA_GLOBAL
void kernel_BaseFabReal(amrex::BaseFab<amrex::Real> *bf1, amrex::Real *val, amrex::IntVect *iv)
{
  amrex::Box bx(bf1->box());
  *iv = bx.bigEnd();
  *val = bf1->nComp();

  bf1->setVal(17.499, bx, 0, 1);
  bf1->getVal(val, amrex::IntVect(7,7,7));
}

AMREX_CUDA_GLOBAL
void kernel_BaseFabInt(amrex::BaseFab<int> *bf1, int *val, amrex::IntVect *iv)
{
  amrex::Box bx(bf1->box());
  *iv = bx.bigEnd();
  *val = bf1->nComp();

  bf1->setVal(17.499, bx, 0, 1);
  bf1->getVal(val, amrex::IntVect(7,7,7));
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
    // BaseFab<Real>
    {
      amrex::Print() << "BaseFab Tests"  << std::endl <<
                        "........................." << std::endl;
 

      amrex::Real val_gpu = 101;
      amrex::Real *d_val_gpu;
      cudaMalloc(&d_val_gpu, size_t(sizeof(amrex::Real)));
      cudaMemcpy(d_val_gpu, &val_gpu, sizeof(amrex::Real), cudaMemcpyHostToDevice);

      int int_gpu = 101;
      int *d_int_gpu;
      cudaMalloc(&d_int_gpu, size_t(sizeof(int)));
      cudaMemcpy(d_int_gpu, &val_gpu, sizeof(int), cudaMemcpyHostToDevice);

      int int_old = -1, int_new = -2;
      amrex::Real val_old = 10, val_new = 11;
      amrex::IntVect *bigEnd = new amrex::IntVect(0,0,0);
      amrex::Box baseBox(amrex::IntVect(0,0,0), amrex::IntVect(14,14,14));

      amrex::BaseFab<amrex::Real> *d_bf1 = new amrex::BaseFab<amrex::Real>(baseBox);
      amrex::BaseFab<int> *d_bf2 = new amrex::BaseFab<int>(baseBox);

      d_bf1->setVal(77.017, baseBox, 0, 1);
      d_bf2->setVal(-77, baseBox, 0, 1);
      d_bf1->getVal(&val_old,amrex::IntVect(1,1,1));
      d_bf2->getVal(&int_old, amrex::IntVect(1,1,1));

      amrex::Print() << "Box = " << baseBox << std::endl;

      cudaPointerAttributes ptr_attr;
      cudaPointerGetAttributes(&ptr_attr, d_bf1);
      amrex::Print() << "d_bf1.isManaged = " << ptr_attr.isManaged << std::endl;
      cudaPointerGetAttributes(&ptr_attr, d_bf1->dataPtr());
      amrex::Print() << "d_bf1.dataPtr().isManaged = " << ptr_attr.isManaged << std::endl;
      cudaPointerGetAttributes(&ptr_attr, d_bf2);
      amrex::Print() << "d_bf2.isManaged = " << ptr_attr.isManaged << std::endl;
      cudaPointerGetAttributes(&ptr_attr, d_bf2->dataPtr());
      amrex::Print() << "d_bf2.dataPtr().isManaged = " << ptr_attr.isManaged << std::endl;

//  amrex::Real value_new = 20.0;
//  amrex_fort_fab_setval(baseBox.loVect(), baseBox.hiVect(), d_bf1->dataPtr(0), baseBox.loVect(), baseBox.hiVect(), 0, &value_new);

      int blockSize = 1;
      int numBlocks = 1;
      kernel_BaseFabReal<<<numBlocks,blockSize>>>(d_bf1, d_val_gpu, bigEnd);
      kernel_BaseFabInt <<<numBlocks,blockSize>>>(d_bf2, d_int_gpu, bigEnd);

      cudaDeviceSynchronize();
      cudaMemcpy(&val_gpu, d_val_gpu, sizeof(amrex::Real), cudaMemcpyDeviceToHost); 
      cudaMemcpy(&int_gpu, d_int_gpu, sizeof(int), cudaMemcpyDeviceToHost); 


      d_bf1->getVal(&val_new,amrex::IntVect(0,0,0));
      d_bf2->getVal(&int_new,amrex::IntVect(1,1,1));

      amrex::Print() << "BaseFab Real Test "  << std::endl <<
                        "........................." << std::endl <<
                        "bigEnd = "      << *bigEnd  << std::endl << 
                        "pre setVal = "  << val_old << std::endl << 
                        "gpu Val = "     << val_gpu << std::endl <<
                        "post setVal = " << val_new << std::endl <<
                        "**********************************" << std::endl << std::endl;

      amrex::Print() << "BaseFab Int Test "  << std::endl <<
                        "........................." << std::endl <<
                        "pre setVal = "  << int_old << std::endl << 
                        "gpu Val = "     << int_gpu << std::endl <<
                        "post setVal = " << int_new << std::endl <<
                        "**********************************" << std::endl << std::endl;

    }

    amrex::Finalize();
}
