#include <set>
#include <random>
#include <limits>
#include <AMReX_Arena.H>
#include <AMReX_BLFort.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>
#include <AMReX_BlockMutex.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_OpenMP.H>

#ifdef AMREX_USE_HIP
#include <hiprand.hpp>
#endif

#ifdef AMREX_USE_HIP
using randState_t = hiprandState_t;
#elif defined(AMREX_USE_CUDA)
using randState_t =  curandState_t;
#elif defined(AMREX_USE_DPCPP)
struct randState_t {};
#endif

namespace
{
    int nthreads;

    amrex::Vector<std::mt19937> generators;

#ifdef AMREX_USE_GPU
    // This seems to be a good default value on NVIDIA V100 GPUs
    constexpr int gpu_nstates_default = 1e5;

    int gpu_nstates_h = 0;
#ifndef AMREX_USE_DPCPP
    AMREX_GPU_DEVICE int gpu_nstates_d = 0;
#endif

    randState_t* d_states_h_ptr = nullptr;
#ifndef AMREX_USE_DPCPP
    AMREX_GPU_DEVICE randState_t* d_states_d_ptr;
#endif

    amrex::BlockMutex* h_mutex_h_ptr = nullptr;
    amrex::BlockMutex* d_mutex_h_ptr = nullptr;

#ifndef AMREX_USE_DPCPP
    AMREX_GPU_DEVICE
    amrex::BlockMutex* d_mutex_d_ptr = nullptr;
#endif

#endif

}

void
amrex::InitRandom (amrex::ULong seed, int nprocs)
{
    nthreads = OpenMP::get_max_threads();
    generators.resize(nthreads);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        int tid = OpenMP::get_thread_num();
        amrex::ULong init_seed = seed + tid*nprocs;
        generators[tid].seed(init_seed);
    }

#ifdef AMREX_USE_GPU
    DeallocateRandomSeedDevArray();
    ResizeRandomSeed(gpu_nstates_default);
#endif
}

#ifdef AMREX_USE_GPU
AMREX_GPU_DEVICE
int amrex::get_state (int tid)
{
#ifdef AMREX_USE_DPCPP
// xxxxx DPCPP todo
    amrex::ignore_unused(tid);
    return 0;
#else
    // block size must evenly divide # of RNG states so we cut off the excess states
    int bsize = blockDim.x * blockDim.y * blockDim.z;
    int nstates = gpu_nstates_d - (gpu_nstates_d % bsize);
    int i = tid % nstates;

    d_mutex_d_ptr->lock(i);

    return i;
#endif
}

AMREX_GPU_DEVICE
void amrex::free_state (int tid)
{
#ifdef AMREX_USE_DPCPP
    amrex::ignore_unused(tid);
// xxxxx DPCPP todo
#else
    int bsize = blockDim.x * blockDim.y * blockDim.z;
    int nstates = gpu_nstates_d - (gpu_nstates_d % bsize);
    int i = tid % nstates;

    d_mutex_d_ptr->unlock(i);
#endif
}
#endif

AMREX_GPU_HOST_DEVICE amrex::Real
amrex::RandomNormal (amrex::Real mean, amrex::Real stddev)
{

    amrex::Real rand;

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;

    int tid = blockId * (blockDim.x * blockDim.y * blockDim.z)
              + (threadIdx.z * (blockDim.x * blockDim.y))
              + (threadIdx.y * blockDim.x) + threadIdx.x ;

    int i = get_state(tid);
#ifdef BL_USE_FLOAT
    AMREX_HIP_OR_CUDA( rand = stddev * hiprand_normal(&d_states_d_ptr[i]) + mean;,
                       rand = stddev *  curand_normal(&d_states_d_ptr[i]) + mean; )
#else
    AMREX_HIP_OR_CUDA( rand = stddev * hiprand_normal_double(&d_states_d_ptr[i]) + mean;,
                       rand = stddev *  curand_normal_double(&d_states_d_ptr[i]) + mean; )
#endif
    __threadfence();
    free_state(tid);

#elif defined(__SYCL_DEVICE_ONLY__)

    amrex::ignore_unused(mean,stddev);
    assert(0);
    rand = Real(0.0);
    return rand;

#else

    std::normal_distribution<amrex::Real> distribution(mean, stddev);
    int tid = OpenMP::get_thread_num();
    rand = distribution(generators[tid]);

#endif
    return rand;
}

AMREX_GPU_HOST_DEVICE amrex::Real
amrex::Random ()
{
    amrex::Real rand;
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)    // on the device
    int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;

    int tid = blockId * (blockDim.x * blockDim.y * blockDim.z)
              + (threadIdx.z * (blockDim.x * blockDim.y))
              + (threadIdx.y * blockDim.x) + threadIdx.x ;
    int i = get_state(tid);

    // curand_uniform generates numbers in (0.0,1], while
    // std::uniform_real_distribution in [0.0, 1.0)
#ifdef BL_USE_FLOAT
    AMREX_HIP_OR_CUDA( rand = 1.0f - hiprand_uniform(&d_states_d_ptr[i]);,
                       rand = 1.0f - curand_uniform(&d_states_d_ptr[i]); )
#else
    AMREX_HIP_OR_CUDA( rand = 1.0 - hiprand_uniform_double(&d_states_d_ptr[i]);,
                       rand = 1.0 - curand_uniform_double(&d_states_d_ptr[i]); )
#endif

    __threadfence();
    free_state(tid);

#elif defined(__SYCL_DEVICE_ONLY__)

    assert(0);
    rand = Real(0.0);
    return rand;

#else     // on the host

    std::uniform_real_distribution<amrex::Real> distribution(0.0, 1.0);
    int tid = OpenMP::get_thread_num();
    rand = distribution(generators[tid]);

#endif

    return rand;
}

AMREX_GPU_HOST_DEVICE unsigned int
amrex::RandomPoisson (amrex::Real lambda)
{
    unsigned int rand;

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    const auto blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;

    const auto tid = blockId * (blockDim.x * blockDim.y * blockDim.z)
              + (threadIdx.z * (blockDim.x * blockDim.y))
              + (threadIdx.y * blockDim.x) + threadIdx.x ;

    const auto i = get_state(tid);

    AMREX_HIP_OR_CUDA( rand = hiprand_poisson(&d_states_d_ptr[i], lambda);,
                       rand = curand_poisson(&d_states_d_ptr[i], lambda);)

    __threadfence();
    free_state(tid);

#elif defined(__SYCL_DEVICE_ONLY__)

    amrex::ignore_unused(lambda);
    assert(0);
    rand = 0;
    return rand;

#else

    std::poisson_distribution<unsigned int> distribution(lambda);
    int tid = OpenMP::get_thread_num();
    rand = distribution(generators[tid]);

#endif
    return rand;
}

AMREX_GPU_HOST_DEVICE unsigned int
amrex::Random_int (unsigned int n)
{
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)  // on the device
    constexpr unsigned int RAND_M = 4294967295; // 2**32-1

    int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;

    int tid = blockId * (blockDim.x * blockDim.y * blockDim.z)
              + (threadIdx.z * (blockDim.x * blockDim.y))
              + (threadIdx.y * blockDim.x) + threadIdx.x ;

    unsigned int rand;
    int i = get_state(tid);
    do {
        AMREX_HIP_OR_CUDA( rand = hiprand(&d_states_d_ptr[i]);,
                           rand =  curand(&d_states_d_ptr[i]); )
    } while (rand > (RAND_M - RAND_M % n));
    __threadfence();
    free_state(tid);

    return rand % n;

#elif defined(__SYCL_DEVICE_ONLY__)

    amrex::ignore_unused(n);
    assert(0);
    return 0;

#else // on the host

    std::uniform_int_distribution<unsigned int> distribution(0, n-1);
    int tid = OpenMP::get_thread_num();
    return distribution(generators[tid]);

#endif
}

AMREX_GPU_HOST amrex::ULong
amrex::Random_long (amrex::ULong n)
{
    std::uniform_int_distribution<amrex::ULong> distribution(0, n-1);
    int tid = OpenMP::get_thread_num();
    return distribution(generators[tid]);
}

void
amrex::SaveRandomState (std::ostream& os)
{
    for (int i = 0; i < nthreads; i++) {
        os << generators[i] << "\n";
    }
}

void
amrex::RestoreRandomState (std::istream& is, int nthreads_old, int nstep_old)
{
    int N = std::min(nthreads, nthreads_old);
    for (int i = 0; i < N; i++)
        is >> generators[i];
    if (nthreads > nthreads_old) {
        const int NProcs = ParallelDescriptor::NProcs();
        const int MyProc = ParallelDescriptor::MyProc();
        for (int i = nthreads_old; i < nthreads; i++) {
            amrex::ULong seed = MyProc+1 + i*NProcs;
            if (std::numeric_limits<amrex::ULong>::max()/(amrex::ULong)(nstep_old+1)
                > static_cast<amrex::ULong>(nthreads*NProcs)) // avoid overflow
            {
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

void amrex::ResetRandomSeed (amrex::ULong seed)
{
    InitRandom(seed);
}

void
amrex::ResizeRandomSeed (int N)
{
    amrex::ignore_unused(N);

    BL_PROFILE("ResizeRandomSeed");

#ifdef AMREX_USE_DPCPP
    // xxxxx DPCPP todo

#elif defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

    if (N <= gpu_nstates_h) return;

    int PrevSize = gpu_nstates_h;
    int SizeDiff = N - PrevSize;

    randState_t* new_data;
    new_data = static_cast<randState_t*>  (The_Device_Arena()->alloc(N*sizeof(randState_t)));

    if (h_mutex_h_ptr != nullptr)
    {
        delete h_mutex_h_ptr;
        h_mutex_h_ptr = nullptr;
    }

    if (d_mutex_h_ptr != nullptr)
    {
        The_Device_Arena()->free(d_mutex_h_ptr);
        d_mutex_h_ptr = nullptr;
    }

    h_mutex_h_ptr = new amrex::BlockMutex(N);
    d_mutex_h_ptr = static_cast<amrex::BlockMutex*> (The_Device_Arena()->alloc(sizeof(amrex::BlockMutex)));
    amrex::Gpu::htod_memcpy(d_mutex_h_ptr, h_mutex_h_ptr, sizeof(amrex::BlockMutex));

    if (d_states_h_ptr != nullptr)
    {
        amrex::Gpu::dtod_memcpy(new_data, d_states_h_ptr, PrevSize*sizeof(randState_t));
        The_Device_Arena()->free(d_states_h_ptr);
    }

    d_states_h_ptr = new_data;
    gpu_nstates_h = N;
    amrex::BlockMutex* d_mutex_h_ptr_local = d_mutex_h_ptr;

    // xxxxx HIP FIX HERE - hipMemcpyToSymbol doesn't work with pointers.
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int)
    {
        d_states_d_ptr = new_data;
        d_mutex_d_ptr = d_mutex_h_ptr_local;
        gpu_nstates_d = N;
    });

    const int MyProc = amrex::ParallelDescriptor::MyProc();
    amrex::ParallelFor (SizeDiff, [=] AMREX_GPU_DEVICE (int idx) noexcept
    {
        amrex::ULong seed = MyProc*1234567ULL + 12345ULL ;
        int seqstart = idx + 10 * idx ;
        int loc = idx + PrevSize;

        AMREX_HIP_OR_CUDA( hiprand_init(seed, seqstart, 0, &d_states_d_ptr[loc]);,
                            curand_init(seed, seqstart, 0, &d_states_d_ptr[loc]); )
    });

#endif

}

void
amrex::DeallocateRandomSeedDevArray ()
{
#ifdef AMREX_USE_GPU
    if (d_states_h_ptr != nullptr)
    {
        The_Device_Arena()->free(d_states_h_ptr);
        d_states_h_ptr = nullptr;
    }

    if (h_mutex_h_ptr != nullptr)
    {
        delete h_mutex_h_ptr;
        h_mutex_h_ptr = nullptr;
    }

    if (d_mutex_h_ptr != nullptr)
    {
        The_Device_Arena()->free(d_mutex_h_ptr);
        d_mutex_h_ptr = nullptr;
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

#if !defined(AMREX_XSDK) && !defined(BL_NO_FORT)
BL_FORT_PROC_DECL(BLUTILINITRAND,blutilinitrand)(const int* sd)
{
    amrex::ULong seed = *sd;
    amrex::InitRandom(seed);
}

BL_FORT_PROC_DECL(BLINITRAND,blinitrand)(const int* sd)
{
    amrex::ULong seed = *sd;
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

    // This is for Fortran, which doesn't have unsigned long.
    amrex::Long amrex_random_int (amrex::Long n)
    {
        return static_cast<amrex::Long>(amrex::Random_int(static_cast<amrex::ULong>(n)));
    }
}
