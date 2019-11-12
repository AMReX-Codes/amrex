#include <set>
#include <random>
#include <AMReX_Arena.H>
#include <AMReX_BLFort.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuDevice.H>

#ifdef AMREX_USE_HIP
#include <hiprand.hpp>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef AMREX_USE_HIP
using randState_t = hiprandState_t;,
#elif defined(AMREX_USE_CUDA)
using randState_t =  curandState_t; 
#endif

namespace
{
    int nthreads;

    amrex::Vector<std::mt19937> generators;

#ifdef AMREX_USE_GPU

    // This seems to be a good default value on NVIDIA V100 GPUs
    constexpr int gpu_nstates_default = 1e5;

    AMREX_GPU_DEVICE int gpu_nstates_d = 0;
    int gpu_nstates_h = 0;

    AMREX_GPU_DEVICE randState_t* states_d_ptr;
    randState_t* states_h_ptr = nullptr;

    AMREX_GPU_DEVICE int* locks_d_ptr;
    int* locks_h_ptr = nullptr;

#endif

}

void
amrex::InitRandom (unsigned long seed, int nprocs)
{

#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#else
    nthreads = 1;
#endif
    generators.resize(nthreads);

#ifdef _OPENMP
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        unsigned long init_seed = seed + tid*nprocs;
        generators[tid].seed(init_seed);
    }
#else
    generators[0].seed(seed);
#endif

#ifdef AMREX_USE_GPU
    DeallocateRandomSeedDevArray();
    ResizeRandomSeed(gpu_nstates_default);
#endif
}

#ifdef AMREX_DEVICE_COMPILE
AMREX_GPU_DEVICE
int amrex::get_state(int tid)
{
  // block size must evenly divide # of RNG states so we cut off the excess states
  int bsize = blockDim.x * blockDim.y * blockDim.z;
  int nstates = gpu_nstates_d - (gpu_nstates_d % bsize);
  int i = tid % nstates;
  if (tid % bsize == 0)
    {
      while (amrex::Gpu::Atomic::CAS(&locks_d_ptr[i],0,1))
        {
          continue;  //Trap 1 thread per block
        }
    }
  __syncthreads(); // Other threads will wait here for master to be freed
  return i;        // This ensures threads move with block - Important for prevolta
}

AMREX_GPU_DEVICE
void amrex::free_state(int tid)
{
  int bsize = blockDim.x * blockDim.y * blockDim.z;
  int nstates = gpu_nstates_d - (gpu_nstates_d % bsize);
  int i = tid % nstates;
  if (tid % bsize == 0)  // we only locked the master thread state. 
  {
        locks_d_ptr[i] = 0;
  }
}
#endif

AMREX_GPU_HOST_DEVICE amrex::Real
amrex::RandomNormal (amrex::Real mean, amrex::Real stddev)
{

    amrex::Real rand;

#ifdef AMREX_DEVICE_COMPILE 
    int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;

    int tid = blockId * (blockDim.x * blockDim.y * blockDim.z)
              + (threadIdx.z * (blockDim.x * blockDim.y))
              + (threadIdx.y * blockDim.x) + threadIdx.x ;

    int i = get_state(tid);
#ifdef BL_USE_FLOAT
    AMREX_HIP_OR_CUDA( rand = stddev * hiprand_normal(&states_d_ptr[i]) + mean;,
                       rand = stddev *  curand_normal(&states_d_ptr[i]) + mean; );
#else
    AMREX_HIP_OR_CUDA( rand = stddev * hiprand_normal_double(&states_d_ptr[i]) + mean;,
                       rand = stddev *  curand_normal_double(&states_d_ptr[i]) + mean; );
#endif
    free_state(tid);
#else

#ifdef _OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif

    std::normal_distribution<amrex::Real> distribution(mean, stddev);
    rand = distribution(generators[tid]);

#endif
    return rand;
}

AMREX_GPU_HOST_DEVICE amrex::Real
amrex::Random ()
{
    amrex::Real rand;
#ifdef AMREX_DEVICE_COMPILE
    int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;

    int tid = blockId * (blockDim.x * blockDim.y * blockDim.z)
              + (threadIdx.z * (blockDim.x * blockDim.y))
              + (threadIdx.y * blockDim.x) + threadIdx.x ;
    int i = get_state(tid);
#ifdef BL_USE_FLOAT
    AMREX_HIP_OR_CUDA( rand = hiprand_uniform(&states_d_ptr[i]);,
                       rand =  curand_uniform(&states_d_ptr[i]); );
#else
    AMREX_HIP_OR_CUDA( rand = hiprand_uniform_double(&states_d_ptr[i]);,
                       rand =  curand_uniform_double(&states_d_ptr[i]); );
#endif
    free_state(tid);

#else

#ifdef _OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif

    std::uniform_real_distribution<amrex::Real> distribution(0.0, 1.0);
    rand = distribution(generators[tid]);
#endif

    return rand;
}

#ifdef AMREX_DEVICE_COMPILE 

AMREX_GPU_DEVICE unsigned int
amrex::Random_int (unsigned int N)
{
    constexpr unsigned int RAND_M = 4294967295; // 2**32-1

    int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;

    int tid = blockId * (blockDim.x * blockDim.y * blockDim.z)
              + (threadIdx.z * (blockDim.x * blockDim.y))
              + (threadIdx.y * blockDim.x) + threadIdx.x ;

    unsigned int rand;
    int i = get_state(tid);
    do {
        AMREX_HIP_OR_CUDA( rand = hiprand(&states_d_ptr[i]);,
                           rand =  curand(&states_d_ptr[i]); );
    } while (rand > (RAND_M - RAND_M % N));
    free_state(tid);

    return rand % N;
}

#else

AMREX_GPU_HOST unsigned long
amrex::Random_int (unsigned long n)
{
#ifdef _OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif
    std::uniform_int_distribution<unsigned long> distribution(0, n-1);
    return distribution(generators[tid]);
}

#endif

void
amrex::SaveRandomState(std::ostream& os)
{
    for (int i = 0; i < nthreads; i++) {
        os << generators[i] << "\n";
    }
}

void
amrex::RestoreRandomState(std::istream& is, int nthreads_old, int nstep_old)
{
    int N = std::min(nthreads, nthreads_old);
    for (int i = 0; i < N; i++)
        is >> generators[i];
    if (nthreads > nthreads_old) {
        const int NProcs = ParallelDescriptor::NProcs();
        const int MyProc = ParallelDescriptor::MyProc();
        for (int i = nthreads_old; i < nthreads; i++) {
            unsigned long seed = MyProc+1 + i*NProcs;
            if (ULONG_MAX/(unsigned long)(nstep_old+1) >static_cast<unsigned long>(nthreads*NProcs)) { // avoid overflow
                seed += nstep_old*nthreads*NProcs;
            }

            generators[i].seed(seed);
        }
    }
}

void
amrex::UniqueRandomSubset (Vector<int> &uSet, int setSize, int poolSize,
                           bool printSet)
{
  if(setSize > poolSize) {
    amrex::Abort("**** Error in UniqueRandomSubset:  setSize > poolSize.");
  }
  std::set<int> copySet;
  Vector<int> uSetTemp;
  while(static_cast<int>(copySet.size()) < setSize) {
    int r(amrex::Random_int(poolSize));
    if(copySet.find(r) == copySet.end()) {
      copySet.insert(r);
      uSetTemp.push_back(r);
    }
  }
  uSet = uSetTemp;
  if(printSet) {
    for(int i(0); i < uSet.size(); ++i) {
        amrex::AllPrint() << "uSet[" << i << "]  = " << uSet[i] << std::endl;
    }
  }
}

void amrex::ResetRandomSeed(unsigned long seed)
{
    InitRandom(seed);
}

void 
amrex::ResizeRandomSeed (int N)
{
    BL_PROFILE("ResizeRandomSeed");

#ifdef AMREX_USE_GPU

    if (N <= gpu_nstates_h) return;

    int PrevSize = gpu_nstates_h;
    int SizeDiff = N - PrevSize;

    randState_t* new_data;
    int* new_mutex;

    new_data  = static_cast<randState_t*> (The_Device_Arena()->alloc(N*sizeof(randState_t)));
    new_mutex = static_cast<int*>         (The_Device_Arena()->alloc(N*sizeof(int)));

    if (states_h_ptr != nullptr) {

        AMREX_ASSERT(locks_h_ptr != nullptr);

        amrex::Gpu::dtod_memcpy(new_data, states_h_ptr, PrevSize*sizeof(randState_t));
        amrex::Gpu::dtod_memcpy(new_mutex, locks_h_ptr, PrevSize*sizeof(int));

        The_Device_Arena()->free(states_h_ptr);
        The_Device_Arena()->free(locks_h_ptr);

    }

    states_h_ptr = new_data;
    locks_h_ptr = new_mutex;
    gpu_nstates_h = N;

    // HIP FIX HERE - hipMemcpyToSymbol doesn't work with pointers. 
    AMREX_GPU_LAUNCH_DEVICE(Gpu::ExecutionConfig(1, 1, 0),
    [=] AMREX_GPU_DEVICE
    {
        states_d_ptr = new_data;
        locks_d_ptr = new_mutex;
        gpu_nstates_d = N;
    });

    const int MyProc = amrex::ParallelDescriptor::MyProc();
    AMREX_PARALLEL_FOR_1D (SizeDiff, idx,
    {
        unsigned long seed = MyProc*1234567UL + 12345UL ;
        int seqstart = idx + 10 * idx ;
        int loc = idx + PrevSize;
        locks_d_ptr[loc] = 0;
        AMREX_HIP_OR_CUDA( hiprand_init(seed, seqstart, 0, &states_d_ptr[loc]);,
                            curand_init(seed, seqstart, 0, &states_d_ptr[loc]); );
    });

#endif

}

void
amrex::DeallocateRandomSeedDevArray()
{
#ifdef AMREX_USE_GPU
    if (states_h_ptr != nullptr)
    {
        The_Device_Arena()->free(states_h_ptr);
        states_h_ptr = nullptr;
    }

    if (locks_h_ptr != nullptr)
    {
        The_Device_Arena()->free(locks_h_ptr);
        locks_h_ptr = nullptr;
    }

    gpu_nstates_h = 0;
#endif
}

void
amrex::NItemsPerBin (int totalItems, Vector<int> &binCounts)
{
  if(binCounts.size() == 0) {
    return;
  }
  bool verbose(false);
  int countForAll(totalItems / binCounts.size());
  int remainder(totalItems % binCounts.size());
  if(verbose) {
      amrex::Print() << "amrex::NItemsPerBin:  countForAll remainder = " << countForAll
                     << "  " << remainder << std::endl;
  }
  for(int i(0); i < binCounts.size(); ++i) {
    binCounts[i] = countForAll;
  }
  for(int i(0); i < remainder; ++i) {
    ++binCounts[i];
  }
  for(int i(0); i < binCounts.size(); ++i) {
    if(verbose) {
        amrex::Print() << "amrex::NItemsPerBin::  binCounts[" << i << "] = " << binCounts[i] << std::endl;
    }
  }
}


//
// Fortran entry points for amrex::Random().
//

#ifndef AMREX_XSDK
BL_FORT_PROC_DECL(BLUTILINITRAND,blutilinitrand)(const int* sd)
{
    unsigned long seed = *sd;
    amrex::InitRandom(seed);
}

BL_FORT_PROC_DECL(BLINITRAND,blinitrand)(const int* sd)
{
    unsigned long seed = *sd;
    amrex::InitRandom(seed);
}

BL_FORT_PROC_DECL(BLUTILRAND,blutilrand)(amrex::Real* rn)
{
    *rn = amrex::Random();
}
#endif

extern "C" {
    double amrex_random ()
    {
        return amrex::Random();
    }

    long amrex_random_int (long n)  // This is for Fortran, which doesn't have unsigned long.
    {
        return static_cast<long>(amrex::Random_int(static_cast<unsigned long>(n)));
    }
}


