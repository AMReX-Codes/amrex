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
#ifdef AMREX_USE_SYCL
    sycl_rng_descr* rand_engine_descr = nullptr;
#else
    amrex::randState_t* gpu_rand_state = nullptr;
#endif
}

namespace {
#ifdef AMREX_USE_SYCL
    oneapi::mkl::rng::philox4x32x10* gpu_rand_generator = nullptr;
#else
    amrex::randGenerator_t gpu_rand_generator = nullptr;
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

#ifdef AMREX_USE_SYCL

    rand_engine_descr = new sycl_rng_descr
        (Gpu::Device::streamQueue(), sycl::range<1>(N), gpu_seed, 1);

    gpu_rand_generator = new std::remove_pointer_t<decltype(gpu_rand_generator)>
        (Gpu::Device::streamQueue(), gpu_seed+1234ULL);

#elif defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

    gpu_rand_state =  static_cast<randState_t*>(The_Arena()->alloc(N*sizeof(randState_t)));
    randState_t* gpu_rand_state_local = gpu_rand_state;
    amrex::ParallelFor(N, [=] AMREX_GPU_DEVICE (int idx) noexcept
    {
        ULong seqstart = static_cast<ULong>(idx) + 10 * static_cast<ULong>(idx);
        AMREX_HIP_OR_CUDA( hiprand_init(gpu_seed, seqstart, 0, &gpu_rand_state_local[idx]);,
                            curand_init(gpu_seed, seqstart, 0, &gpu_rand_state_local[idx]); )
    });

#if defined(AMREX_USE_CUDA)
    AMREX_CURAND_SAFE_CALL(curandCreateGenerator
                           (&gpu_rand_generator, CURAND_RNG_PSEUDO_DEFAULT));
    AMREX_CURAND_SAFE_CALL(curandSetPseudoRandomGeneratorSeed
                           (gpu_rand_generator, gpu_seed+1234ULL));
#else
    AMREX_HIPRAND_SAFE_CALL(hiprandCreateGenerator
                            (&gpu_rand_generator, HIPRAND_RNG_PSEUDO_DEFAULT));
    AMREX_HIPRAND_SAFE_CALL(hiprandSetPseudoRandomGeneratorSeed
                            (gpu_rand_generator, gpu_seed+1234ULL));
#endif

#endif

    Gpu::synchronize();
}
}
#endif

namespace amrex {

void
InitRandom (ULong cpu_seed, int nprocs, ULong gpu_seed)
{
    nthreads = OpenMP::get_max_threads();
    generators.resize(nthreads);

#ifdef AMREX_USE_OMP
    if (omp_in_parallel()) {
        amrex::Abort("It is not safe to call amrex::InitRandom inside a threaded region.");
    }
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    {
        int tid = OpenMP::get_thread_num();
        ULong init_seed = cpu_seed + tid*nprocs;
        generators[tid].seed(init_seed);
    }

#ifdef AMREX_USE_GPU
    ResizeRandomSeed(gpu_seed);
#else
    ignore_unused(gpu_seed);
#endif
}

Real RandomNormal (Real mean, Real stddev)
{
    std::normal_distribution<Real> distribution(mean, stddev);
    int tid = OpenMP::get_thread_num();
    return distribution(generators[tid]);
}

Real Random ()
{
    std::uniform_real_distribution<Real> distribution(0.0, 1.0);
    int tid = OpenMP::get_thread_num();
    return distribution(generators[tid]);
}

unsigned int RandomPoisson (Real lambda)
{
    std::poisson_distribution<unsigned int> distribution(lambda);
    int tid = OpenMP::get_thread_num();
    return distribution(generators[tid]);
}

Real RandomGamma (Real alpha, Real beta)
{
    std::gamma_distribution<Real> distribution(alpha, beta);
    int tid = OpenMP::get_thread_num();
    return distribution(generators[tid]);
}

unsigned int Random_int (unsigned int n)
{
    std::uniform_int_distribution<unsigned int> distribution(0, n-1);
    int tid = OpenMP::get_thread_num();
    return distribution(generators[tid]);
}

ULong Random_long (ULong n)
{
    std::uniform_int_distribution<ULong> distribution(0, n-1);
    int tid = OpenMP::get_thread_num();
    return distribution(generators[tid]);
}

void
SaveRandomState (std::ostream& os)
{
    for (int i = 0; i < nthreads; i++) {
        os << generators[i] << "\n";
    }
}

void
RestoreRandomState (std::istream& is, int nthreads_old, int nstep_old)
{
    int N = std::min(nthreads, nthreads_old);
    for (int i = 0; i < N; i++) {
        is >> generators[i];
    }
    if (nthreads > nthreads_old) {
        const int NProcs = ParallelDescriptor::NProcs();
        const int MyProc = ParallelDescriptor::MyProc();
        for (int i = nthreads_old; i < nthreads; i++) {
            ULong seed = MyProc+1 + i*NProcs;
            if (std::numeric_limits<ULong>::max()/static_cast<ULong>(nstep_old+1)
                > static_cast<ULong>(nthreads)*static_cast<ULong>(NProcs)) // avoid overflow
            {
                seed += nstep_old*nthreads*NProcs;
            }

            generators[i].seed(seed);
        }
    }
}

void
UniqueRandomSubset (Vector<int> &uSet, int setSize, int poolSize,
                    bool printSet)
{
  if(setSize > poolSize) {
    Abort("**** Error in UniqueRandomSubset:  setSize > poolSize.");
  }
  std::set<int> copySet;
  uSet.clear();
  while(static_cast<int>(copySet.size()) < setSize) {
    int r = static_cast<int>(Random_int(poolSize));
    if(copySet.find(r) == copySet.end()) {
      copySet.insert(r);
      uSet.push_back(r);
    }
  }
  if(printSet) {
    for(int i(0); i < uSet.size(); ++i) {
        AllPrint() << "uSet[" << i << "]  = " << uSet[i] << '\n';
    }
  }
}

void ResetRandomSeed (ULong cpu_seed, ULong gpu_seed)
{
    InitRandom(cpu_seed, ParallelDescriptor::NProcs(), gpu_seed);
}

void
DeallocateRandomSeedDevArray ()
{
#ifdef AMREX_USE_GPU
#ifdef AMREX_USE_SYCL
    if (rand_engine_descr) {
        delete rand_engine_descr;
        Gpu::streamSynchronize();
        rand_engine_descr = nullptr;
    }
    if (gpu_rand_generator != nullptr) {
        delete gpu_rand_generator;
        Gpu::streamSynchronize();
        gpu_rand_generator = nullptr;
    }
#else
    if (gpu_rand_state != nullptr)
    {
        The_Arena()->free(gpu_rand_state);
        gpu_rand_state = nullptr;
    }
    if (gpu_rand_generator != nullptr)
    {
#if defined(AMREX_USE_CUDA)
        AMREX_CURAND_SAFE_CALL(curandDestroyGenerator(gpu_rand_generator));
#else
        AMREX_HIPRAND_SAFE_CALL(hiprandDestroyGenerator(gpu_rand_generator));
#endif
        gpu_rand_generator = nullptr;
    }
#endif
#endif
}

void FillRandom (Real* p, Long N)
{
#ifdef AMREX_USE_CUDA

#  ifdef BL_USE_FLOAT
    AMREX_CURAND_SAFE_CALL(curandGenerateUniform(gpu_rand_generator, p, N));
#  else
    AMREX_CURAND_SAFE_CALL(curandGenerateUniformDouble(gpu_rand_generator, p, N));
#  endif
    Gpu::synchronize();

#elif defined(AMREX_USE_HIP)

#  ifdef BL_USE_FLOAT
    AMREX_HIPRAND_SAFE_CALL(hiprandGenerateUniform(gpu_rand_generator, p, N));
#  else
    AMREX_HIPRAND_SAFE_CALL(hiprandGenerateUniformDouble(gpu_rand_generator, p, N));
#  endif
    Gpu::synchronize();

#elif defined(AMREX_USE_SYCL)

    oneapi::mkl::rng::uniform<Real> distr;
    auto event = oneapi::mkl::rng::generate(distr, *gpu_rand_generator, N, p);
    event.wait();

#else
    std::uniform_real_distribution<Real> distribution(Real(0.0), Real(1.0));
    auto& gen = generators[OpenMP::get_thread_num()];
    for (Long i = 0; i < N; ++i) {
        p[i] = distribution(gen);
    }
#endif
}

void FillRandomNormal (Real* p, Long N, Real mean, Real stddev)
{
    if (N <= 0) { return; }

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    if (N == 1) {
        auto r = amrex::RandomNormal(mean, stddev);
        Gpu::htod_memcpy_async(p, &r, sizeof(Real));
        Gpu::streamSynchronize();
        return;
    }
    // The length passed to [cu|hip]randGenerateNormal must be even
    Long Neven =  (N%2 == 0) ? N : N-1;
#endif

#if defined(AMREX_USE_CUDA)

#  ifdef BL_USE_FLOAT
    AMREX_CURAND_SAFE_CALL(curandGenerateNormal(gpu_rand_generator, p, Neven, mean, stddev));
#  else
    AMREX_CURAND_SAFE_CALL(curandGenerateNormalDouble(gpu_rand_generator, p, Neven, mean, stddev));
#  endif

#elif defined(AMREX_USE_HIP)

#  ifdef BL_USE_FLOAT
    AMREX_HIPRAND_SAFE_CALL(hiprandGenerateNormal(gpu_rand_generator, p, Neven, mean, stddev));
#  else
    AMREX_HIPRAND_SAFE_CALL(hiprandGenerateNormalDouble(gpu_rand_generator, p, Neven, mean, stddev));
#  endif

#elif defined(AMREX_USE_SYCL)

    oneapi::mkl::rng::gaussian<Real> distr(mean, stddev);
    auto event = oneapi::mkl::rng::generate(distr, *gpu_rand_generator, N, p);
    event.wait();

#else

    std::normal_distribution<Real> distribution(mean, stddev);
    auto& gen = generators[OpenMP::get_thread_num()];
    for (Long i = 0; i < N; ++i) {
        p[i] = distribution(gen);
    }

#endif

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    if (Neven < N) {
        auto r = amrex::RandomNormal(mean, stddev);
        Gpu::htod_memcpy_async(p+(N-1), &r, sizeof(Real));
    }
    Gpu::synchronize();
#endif
}

} // namespace amrex

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
