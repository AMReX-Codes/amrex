#include <AMReX_Arena.H>
#include <AMReX_BLFort.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>
#include <AMReX_Gpu.H>
#include <AMReX_OpenMP.H>

#include <set>
#include <random>
#include <limits>

namespace
{
    int nthreads;
    amrex::Vector<std::mt19937> generators;
}

#ifdef AMREX_USE_GPU
namespace amrex {
#ifdef AMREX_USE_DPCPP
    dpcpp_rng_descr* rand_engine_descr = nullptr;
#else
    amrex::randState_t* gpu_rand_state = nullptr;
#endif
}
#endif

#ifdef AMREX_USE_GPU
namespace {
void ResizeRandomSeed (amrex::ULong gpu_seed)
{
    BL_PROFILE("ResizeRandomSeed");

    using namespace amrex;

    DeallocateRandomSeedDevArray();

    const int N = Gpu::Device::maxBlocksPerLaunch() * AMREX_GPU_MAX_THREADS;

#ifdef AMREX_USE_DPCPP

    rand_engine_descr = new dpcpp_rng_descr
        (Gpu::Device::streamQueue(), sycl::range<1>(N), gpu_seed, 1);

#elif defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

    gpu_rand_state =  static_cast<randState_t*>(The_Arena()->alloc(N*sizeof(randState_t)));
    randState_t* gpu_rand_state_local = gpu_rand_state;
    amrex::ParallelFor(N, [=] AMREX_GPU_DEVICE (int idx) noexcept
    {
        ULong seqstart = static_cast<ULong>(idx) + 10 * static_cast<ULong>(idx);
        AMREX_HIP_OR_CUDA( hiprand_init(gpu_seed, seqstart, 0, &gpu_rand_state_local[idx]);,
                            curand_init(gpu_seed, seqstart, 0, &gpu_rand_state_local[idx]); )
    });
#endif

    Gpu::streamSynchronize();
}
}
#endif

void
amrex::InitRandom (amrex::ULong cpu_seed, int nprocs, amrex::ULong gpu_seed)
{
    nthreads = OpenMP::get_max_threads();
    generators.resize(nthreads);

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    {
        int tid = OpenMP::get_thread_num();
        amrex::ULong init_seed = cpu_seed + tid*nprocs;
        generators[tid].seed(init_seed);
    }

#ifdef AMREX_USE_GPU
    ResizeRandomSeed(gpu_seed);
#else
    amrex::ignore_unused(gpu_seed);
#endif
}

amrex::Real amrex::RandomNormal (amrex::Real mean, amrex::Real stddev)
{
    std::normal_distribution<amrex::Real> distribution(mean, stddev);
    int tid = OpenMP::get_thread_num();
    return distribution(generators[tid]);
}

amrex::Real amrex::Random ()
{
    std::uniform_real_distribution<amrex::Real> distribution(0.0, 1.0);
    int tid = OpenMP::get_thread_num();
    return distribution(generators[tid]);
}

unsigned int amrex::RandomPoisson (amrex::Real lambda)
{
    std::poisson_distribution<unsigned int> distribution(lambda);
    int tid = OpenMP::get_thread_num();
    return distribution(generators[tid]);
}

unsigned int amrex::Random_int (unsigned int n)
{
    std::uniform_int_distribution<unsigned int> distribution(0, n-1);
    int tid = OpenMP::get_thread_num();
    return distribution(generators[tid]);
}

amrex::ULong amrex::Random_long (amrex::ULong n)
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

void amrex::ResetRandomSeed (amrex::ULong cpu_seed, amrex::ULong gpu_seed)
{
    InitRandom(cpu_seed, ParallelDescriptor::NProcs(), gpu_seed);
}

void
amrex::DeallocateRandomSeedDevArray ()
{
#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_DPCPP
    if (rand_engine_descr) {
        delete rand_engine_descr;
        Gpu::streamSynchronize();
        rand_engine_descr = nullptr;
    }
#else
    if (gpu_rand_state != nullptr)
    {
        The_Arena()->free(gpu_rand_state);
        gpu_rand_state = nullptr;
    }
#endif
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
