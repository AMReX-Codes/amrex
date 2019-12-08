#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <AMReX_Utility.H>
#include <AMReX_Array.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BlockMutex.H>

using namespace amrex;

void lockingTest(int num_draw);
void blockCountingTest();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        int num_draw;
        
        ParmParse pp;        
        pp.get("num_draw", num_draw);
                    
        auto begin = std::chrono::high_resolution_clock::now();
        amrex::Print() << "Testing using locks to do atomic adds \n";
        for (int i = 0; i < 10; ++i)
            lockingTest(num_draw);
        amrex::Print() << "Testing using locks to count the number of blocks a kernel was launched with \n";
        for (int i = 0; i < 10; ++i)
            blockCountingTest();
        amrex::Print() << "Locking test passed! \n";
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Execution Time: ";
        std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()
                  << "ns" << std::endl;
    }
    amrex::Finalize();
}

AMREX_GPU_DEVICE
void addone(double volatile * num)
{
  int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;
  
  int tid = blockId * (blockDim.x * blockDim.y * blockDim.z)
    + (threadIdx.z * (blockDim.x * blockDim.y))
    + (threadIdx.y * blockDim.x) + threadIdx.x ;

  int i = amrex::get_state(tid);
  num[i] = num[i]+1;
  __threadfence();
  amrex::free_state(tid);
}

void lockingTest (int Ndraw)
{
#ifdef AMREX_USE_CUDA    
    Gpu::DeviceVector<double> d_xpos(Ndraw, 0.0);
    Gpu::DeviceVector<double> d_ypos(Ndraw, 0.0);
    Gpu::DeviceVector<double> d_zpos(Ndraw, 0.0);

    double *dxpos = d_xpos.dataPtr();
    double *dypos = d_ypos.dataPtr();
    double *dzpos = d_zpos.dataPtr();    
#endif

    amrex::Vector<double> hx(Ndraw);
    amrex::Vector<double> hy(Ndraw);
    amrex::Vector<double> hz(Ndraw);
    double *hxpos = hx.dataPtr();
    double *hypos = hy.dataPtr();
    double *hzpos = hz.dataPtr();

    AMREX_PARALLEL_FOR_1D (Ndraw, i,
    {
#ifdef AMREX_USE_CUDA
        addone(dxpos);
        addone(dypos);
        addone(dzpos);	  
#else // Not defined for CPU -> this will still have the test pass if run with cpus
        hx[i] = hx[i]+1;
        hy[i] = hy[i]+1;
        hz[i] = hz[i]+1;
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

   // For the test to pass and race conditions avoided successfuly
   // All of these numbers should be the same
   amrex::Print() << "Sumx  = " << sumx << std::endl;
   amrex::Print() << "Sumy  = " << sumy << std::endl;
   amrex::Print() << "Sumz  = " << sumz << std::endl;
   amrex::Print() << "Ndraw = " << Ndraw << std::endl;
   amrex::Print() << "\n";
   
   AMREX_ALWAYS_ASSERT((sumx == sumy) && (sumz == Ndraw) && (sumx == sumz));   
}

__global__
void count_blocks (BlockMutex* mut, int *numBlocks, int lock_index)
{
    mut->lock(lock_index);
    if (threadIdx.x == 0) { numBlocks[0] = numBlocks[0] + 1; }
    __threadfence();
    mut->unlock(lock_index);
}

void blockCountingTest ()
{
    
    constexpr int NUMBLOCKS = 512;
    constexpr int NUMTHREADS = 1024;
    
    {
        int h_counting, *d_counting;
        
        AMREX_CUDA_SAFE_CALL(cudaMalloc(&d_counting, sizeof(int)));

        h_counting = 0;

        AMREX_CUDA_SAFE_CALL(cudaMemcpy(d_counting, &h_counting, sizeof(int), cudaMemcpyHostToDevice));

        BlockMutex h_mut(1);
        BlockMutex* d_mut;
        AMREX_CUDA_SAFE_CALL(cudaMalloc(&d_mut, sizeof(BlockMutex)));
        AMREX_CUDA_SAFE_CALL(cudaMemcpy(d_mut, &h_mut, sizeof(BlockMutex), cudaMemcpyHostToDevice));
               
        count_blocks<<<NUMBLOCKS, NUMTHREADS>>>(d_mut, d_counting, 0);

        cudaPeekAtLastError();

        Gpu::Device::synchronize();

        AMREX_CUDA_SAFE_CALL(cudaMemcpy(&h_counting, d_counting, sizeof(int), cudaMemcpyDeviceToHost));

        amrex::Print() << "Number of blocks is: " << h_counting << "\n";

        if (h_counting != NUMBLOCKS) amrex::Abort();

        AMREX_CUDA_SAFE_CALL(cudaFree(d_counting));
    }
}
