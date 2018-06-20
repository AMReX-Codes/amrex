#include <cuda_device_runtime_api.h>

#include <iostream>
#include <AMReX.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_BaseFab.H>
#include <AMReX_BaseFab_f.H>

#define RLIM_3D(x) &((int []){x[0], x[1], 0, x[2], 0}[0])

AMREX_CUDA_DEVICE
void myf(int* a)
{
   printf("a[1] = %i", a[1]);
}

AMREX_CUDA_GLOBAL
void passBoxByValue(amrex::Box bx)
{
   printf("Box.bigEnd[2] = %i\n", bx.bigEnd(2));

   myf(RLIM_3D(bx.bigEnd()));

//   printf("array[0] = %i", arr[2]);
}


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

   amrex::BaseFab<amrex::Real> local(bx);
   local.setVal(43.001, bx, 0, 1);
   local.getVal(val, amrex::IntVect(7,7,7));
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

AMREX_CUDA_GLOBAL
void kernel_BaseFabCopy(amrex::BaseFab<amrex::Real> *bf1, amrex::BaseFab<amrex::Real> *bf2)
{
   amrex::Real* mem = (new amrex::Real(1000));
//  cudaMalloc(&mem, 1000);
   amrex::Box bx1 = bf1->box(); 
   amrex::Box bx2 = bx1;
   for (int i=0; i<3; ++i)
   {
     bx2.setBig(i,(bx2.bigEnd(i))/2);
   }

   bf1->copy(*bf2, bx1, 0, bx2, 0, 1);

   bf1->copyToMem(bx1, 0, 1, mem);
   bf2->copyFromMem(bx1, 0, 1, mem);
   delete(mem);
//  cudaFree(mem);
}
/*
AMREX_CUDA_GLOBAL
void kernel_Geometry(amrex::GeometryData *geom, amrex::Box *domain, amrex::Real *off, amrex::Box *bx)
{
    for (int i=0; i<3; ++i)
    {
      off[i] = geom->offset[i];
    }

    *bx = *domain;
}
*/
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

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

      amrex::BaseFab<amrex::Real> *d_bf3 = new amrex::BaseFab<amrex::Real>(baseBox);
      d_bf3->setVal(-0.017, baseBox, 0, 1);
     

      kernel_BaseFabCopy<<<numBlocks,blockSize>>>(d_bf1, d_bf3);
      cudaDeviceSynchronize();
      d_bf1->getVal(&val_old,d_bf1->smallEnd());
      d_bf1->getVal(&val_new,d_bf1->bigEnd());

      amrex::Print() << "BaseFab Copy Test " << std::endl <<
                        "........................." << std::endl <<
                        "bf1 Min Val = "  << val_old << std::endl << 
                        "bf1 Max Val = "  << val_new << std::endl; 

      d_bf3->getVal(&val_old,d_bf3->smallEnd());
      d_bf3->getVal(&val_new,d_bf3->bigEnd());

      amrex::Print() << "bf2 Min Val = "  << val_old << std::endl << 
                        "bf2 Max Val = "  << val_new << std::endl; 

    }
/*
    // GeomData test
    {
      amrex::Print() << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl << std::endl;

      amrex::Geometry geom;
      amrex::Box baseBox(amrex::IntVect(0,0,0), amrex::IntVect(14,14,14));
      amrex::Box *baseBoxB = new amrex::Box(amrex::IntVect(9,9,9), amrex::IntVect(34,34,34));
      amrex::RealBox *real = new amrex::RealBox(-10.0, -10.0, 0, 25, 25, 35);
      amrex::RealBox *realB = new amrex::RealBox(0, 0, 0, 0, 0, 0);
      int coordsys = 0;    // 0 (1D, 2D & 3D), 1(2D), 2(1D)
      int is_per[3] = {0,0,0};

      geom.define(baseBox, real, coordsys, is_per); 

      amrex::Real *off;
      cudaMallocManaged(&off, size_t(3*sizeof(amrex::Real)));
      amrex::GeometryData *gData = geom.dataPtr();
      for (int i=0; i<3; ++i)
      {
        off[i] = i;
        amrex::Print() << "off[" << i << "] = " << off[i] << std::endl;
        amrex::Print() << "offset[" << i << "] = " << gData->offset[i] << std::endl;
      } 

      amrex::Print() << "off.isManaged = " << checkManaged(off) << std::endl;
      amrex::Print() << "geom.dataPtr().isManaged = " << checkManaged(geom.dataPtr()) << std::endl;
      amrex::Print() << "&(geom.dataPtr()->offset).isManaged = " << checkManaged(&(gData->offset[0])) << std::endl;
      amrex::Print() << "geom.dataPtr()->prob_domain.lo() = " << checkManaged(geom.dataPtr()->prob_domain.lo()) << std::endl;
      amrex::Print() << "geom.dataPtr()->domain = " << gData->domain << std::endl;

      amrex::Print() << "geom.dataPtr()->c_sys = " << gData->c_sys << std::endl; 

      amrex::Print() << "RealBox Before = " << *baseBoxB << std::endl;
      kernel_Geometry<<<1, 1>>>(gData, &(geom.dataPtr()->domain), off, baseBoxB);
      cudaDeviceSynchronize();
      amrex::Print() << "RealBox After = " << *baseBoxB << std::endl;

      for (int i=0; i<3; ++i)
      {
        amrex::Print() << "off[" << i << "] = " << off[i] << std::endl;
      } 

      cudaFree(off);
    }
*/
    {
      amrex::Box vbx(amrex::IntVect(-12,-13,-14), amrex::IntVect(23,24,25));
      passBoxByValue<<<1,1>>>(vbx);
      cudaDeviceSynchronize();
    }

    amrex::Print() << std::endl << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << std::endl << std::endl;

    amrex::Finalize();
}
