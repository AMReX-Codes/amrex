#include <set>
#include <random>
#include <AMReX_BLFort.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
    int nthreads;

    amrex::Vector<std::mt19937> generators;

#ifdef AMREX_USE_CUDA

    // This seems to be a good default value on NVIDIA V100 GPUs
    constexpr int cuda_nstates_default = 1e5;

    AMREX_GPU_DEVICE_MANAGED int cuda_nstates = 0;

    AMREX_GPU_DEVICE curandState_t* states_d_ptr;
    curandState_t* states_h_ptr = nullptr;

    AMREX_GPU_DEVICE int* locks_d_ptr;
    int* locks_h_ptr = nullptr;

#endif

}

// HIP FIX HERE - REWRITE RANDOM FOR HIP AS WELL (hiprand)
// https://github.com/ROCm-Developer-Tools/HIP/blob/master/docs/markdown/CURAND_API_supported_by_HIP.md

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

#ifdef AMREX_USE_CUDA
    DeallocateRandomSeedDevArray();
    ResizeRandomSeed(cuda_nstates_default);
#endif
}

#ifdef __CUDA_ARCH__
AMREX_GPU_DEVICE
int amrex::get_state(int tid)
{
  // block size must evenly divide # of RNG states so we cut off the excess states
  int bsize = blockDim.x * blockDim.y * blockDim.z;
  int nstates = cuda_nstates - (cuda_nstates % bsize);
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
  int nstates = cuda_nstates - (cuda_nstates % bsize);
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

#ifdef __CUDA_ARCH__
    int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;

    int tid = blockId * (blockDim.x * blockDim.y * blockDim.z)
              + (threadIdx.z * (blockDim.x * blockDim.y))
              + (threadIdx.y * blockDim.x) + threadIdx.x ;

    int i = get_state(tid);
#ifdef BL_USE_FLOAT
    rand = stddev * curand_normal(&states_d_ptr[i]) + mean;
#else
    rand = stddev * curand_normal_double(&states_d_ptr[i]) + mean;
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
#ifdef __CUDA_ARCH__
    int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;

    int tid = blockId * (blockDim.x * blockDim.y * blockDim.z)
              + (threadIdx.z * (blockDim.x * blockDim.y))
              + (threadIdx.y * blockDim.x) + threadIdx.x ;
    int i = get_state(tid);
#ifdef BL_USE_FLOAT
    rand = curand_uniform(&states_d_ptr[i]);
#else
    rand = curand_uniform_double(&states_d_ptr[i]);
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

#ifdef __CUDA_ARCH__

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
        rand = curand(&states_d_ptr[i]);
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

#ifdef AMREX_USE_CUDA  

    if (N <= cuda_nstates) return;

    int PrevSize = cuda_nstates;
    int SizeDiff = N - PrevSize;

    curandState_t* new_data;
    int* new_mutex;

    AMREX_CUDA_SAFE_CALL(cudaMalloc(&new_data, N*sizeof(curandState_t)));
    AMREX_CUDA_SAFE_CALL(cudaMalloc(&new_mutex, N*sizeof(int)));

    if (states_h_ptr != nullptr) {

        AMREX_ASSERT(locks_h_ptr != nullptr);

        AMREX_CUDA_SAFE_CALL(cudaMemcpy(new_data, states_h_ptr,
                                        PrevSize*sizeof(curandState_t),
                                        cudaMemcpyDeviceToDevice));
        AMREX_CUDA_SAFE_CALL(cudaMemcpy(new_mutex, locks_h_ptr,
                                        PrevSize*sizeof(int),
                                        cudaMemcpyDeviceToDevice));

        AMREX_CUDA_SAFE_CALL(cudaFree(states_h_ptr));
        AMREX_CUDA_SAFE_CALL(cudaFree(locks_h_ptr));
    }

    states_h_ptr = new_data;
    locks_h_ptr = new_mutex;
    cuda_nstates = N;

    AMREX_CUDA_SAFE_CALL(cudaMemcpyToSymbol(states_d_ptr,
                                            &states_h_ptr,
                                            sizeof(curandState_t *)));
    AMREX_CUDA_SAFE_CALL(cudaMemcpyToSymbol(locks_d_ptr,
                                            &locks_h_ptr,
                                            sizeof(int *)));

    const int MyProc = amrex::ParallelDescriptor::MyProc();
    AMREX_PARALLEL_FOR_1D (SizeDiff, idx,
    {
        unsigned long seed = MyProc*1234567UL + 12345UL ;
        int seqstart = idx + 10 * idx ;
        int loc = idx + PrevSize;
        locks_d_ptr[loc] = 0;
        curand_init(seed, seqstart, 0, &states_d_ptr[loc]);
    });

#endif
}

void
amrex::DeallocateRandomSeedDevArray()
{
#ifdef AMREX_USE_CUDA  
    if (states_h_ptr != nullptr)
    {
        cudaFree(states_h_ptr);
        states_h_ptr = nullptr;
    }

    if (locks_h_ptr != nullptr)
    {
        cudaFree(locks_h_ptr);
        locks_h_ptr = nullptr;
    }
    cuda_nstates = 0;
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


