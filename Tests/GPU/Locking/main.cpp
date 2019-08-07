#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <AMReX_Utility.H>
#include <AMReX_Array.H>
#include <AMReX_CudaContainers.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void LockingTest();

int main (int argc, char* argv[])
{

    amrex::Initialize(argc,argv);
    LockingTest();
    amrex::Finalize();

}

void LockingTest ()
{
    int Ndraw=1e7;

#ifdef AMREX_USE_CUDA    
    
    amrex::Print() << "Generating random numbers using GPU ";
    amrex::Print() << amrex::Gpu::Device::deviceId() << " on rank ";
    amrex::Print() << amrex::ParallelDescriptor::MyProc() << "\n";

    Gpu::DeviceVector<double> d_xpos(Ndraw);
    Gpu::DeviceVector<double> d_ypos(Ndraw);
    Gpu::DeviceVector<double> d_zpos(Ndraw);

    double *dxpos = d_xpos.dataPtr();
    double *dypos = d_ypos.dataPtr();
    double *dzpos = d_zpos.dataPtr();
    
#else
    amrex::InitRandom(1024UL,1);
#endif

    amrex::Vector<double> hx(Ndraw);
    amrex::Vector<double> hy(Ndraw);
    amrex::Vector<double> hz(Ndraw);
    double *hxpos = hx.dataPtr();
    double *hypos = hy.dataPtr();
    double *hzpos = hz.dataPtr();



    AMREX_PARALLEL_FOR_1D (Ndraw, idx,
	 {
#ifdef AMREX_USE_CUDA
	   dxpos[idx] = dxpos[idx] + 1;
	   dypos[idx] = dypos[idx] + 1;
	   dzpos[idx] = dzpos[idx] + 1;
#else
	   hx[idx] = hx[idx]+1;
	   hy[idx] = hy[idx]+1;
	   hz[idx] = hz[idx]+1;
#endif	   
	 });

    
#ifdef AMREX_USE_CUDA
        cudaMemcpy(hxpos,dxpos,sizeof(double)*Ndraw,cudaMemcpyDeviceToHost);
        cudaMemcpy(hypos,dypos,sizeof(double)*Ndraw,cudaMemcpyDeviceToHost);
        cudaMemcpy(hzpos,dzpos,sizeof(double)*Ndraw,cudaMemcpyDeviceToHost);
#endif

   int sumx = 0;
   int sumy = 0;
   int sumz = 0;
   for (int i = 0; i < Ndraw; i++ )
     {
       sumx = sumx + hxpos[i];
       sumy = sumy + hypos[i];
       sumz = sumz + hzpos[i];
     }
   //For the test to pass and race conditions avoided successfuly
   // All of these numbers should be the same
   amrex::Print() << "Sumx  = " << sumx << std::endl;
   amrex::Print() << "Sumy  = " << sumy << std::endl;
   amrex::Print() << "Sumz  = " << sumz << std::endl;
   amrex::Print() << "Ndraw = " << Ndraw << std::endl;

   if ( (sumx == sumy) && (sumz == Ndraw) && (sumx == sumz) )
     {
       amrex::Print() << "Locking Test Passed!" << std::endl;
     }
   else
     {
       amrex::Print() << "Locking Test Failed!" << std::endl;
     }
  
}





