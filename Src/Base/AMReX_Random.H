#ifndef AMREX_RAND_H
#define AMREX_RAND_H
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_RandomEngine.H>
#include <limits>
#include <cmath>
#include <cstdint>

namespace amrex
{
    /// \cond DOXYGEN_IGNORE
    namespace random_util {

    //! Convert a (0,1] draw to the documented half-open [0,1) by relocating
    //! the single endpoint that is out of range.  This is the whole interval
    //! conversion: no arithmetic is involved, so nothing can be reassociated,
    //! contracted or flushed, and no resolution is lost.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real one_to_zero (Real v) noexcept { return (v == Real(1)) ? Real(0) : v; }

    //! Clamp a value into the half-open unit interval [0,1) by pulling an
    //! exact 1.0 down to the largest representable Real below one.  Used
    //! where the generator is natively [0,1) and only needs guarding against
    //! reaching its own upper bound; unlike one_to_zero this does not raise
    //! the probability of returning zero.  A comparison rather than
    //! arithmetic, so it holds under every floating-point mode.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real clamp_below_one (Real v) noexcept
    {
        // The constexpr equivalent of std::nextafter(Real(1), Real(0)), which
        // is not available in device code.
        constexpr Real almost_one
            = Real(1) - Real(0.5)*std::numeric_limits<Real>::epsilon();
        // No generator produces NaN, but note this maps it to almost_one
        // rather than propagating it, since the comparison is false for NaN.
        return (v < almost_one) ? v : almost_one;
    }

    }
    /// \endcond

    /**
    * \ingroup amrex_utilities
    * \brief Generate a psuedo-random real from uniform distribution
    *
    *  Generates one pseudorandom real number from a uniform
    *  distribution between 0.0 and 1.0 (0.0 included, 1.0 excluded).
    *
    *  \note Since 0.0 is included, the result must not be passed unguarded
    *  to a function that is singular at zero, such as std::log or std::pow
    *  with a negative exponent.  Note that `1 - amrex::Random()` is \b not
    *  a way around that: it costs an arithmetic operation and discards
    *  resolution, because the subtraction is exact only above 0.5 and
    *  collapses smaller values onto a coarse grid.  Before the upper bound
    *  above was guaranteed it was unsafe as well, evaluating to exactly 0.0
    *  whenever Random() returned exactly 1.0.
    *
    */
    Real Random ();

    /// \ingroup amrex_utilities
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real Random (RandomEngine const& random_engine)
    {
#if defined (__SYCL_DEVICE_ONLY__)
        // clamped for the same reason as the host path: the documented [0,1)
        // of the vendor generator is not a bound we want to depend on
        mkl::rng::device::uniform<Real> distr;
        return random_util::clamp_below_one(mkl::rng::device::generate(distr, *random_engine.engine));
#else
        // hiprand/curand draw from (0,1].  Relocating the single endpoint
        // that is out of range converts that to the documented [0,1) exactly,
        // with no arithmetic.  Subtracting from one would also convert the
        // interval, but not exactly: for a draw below half an ULP of one,
        // `1 - draw` rounds back up to 1.0, and the subtraction additionally
        // collapses every value below 0.5 onto a 2^-24 grid.
#ifdef BL_USE_FLOAT
        AMREX_IF_ON_DEVICE((
                AMREX_HIP_OR_CUDA(
                        return random_util::one_to_zero(hiprand_uniform(random_engine.rand_state)); ,
                        return random_util::one_to_zero(curand_uniform(random_engine.rand_state));
                )
        ))
#else
        AMREX_IF_ON_DEVICE((
                AMREX_HIP_OR_CUDA(
                        return random_util::one_to_zero(hiprand_uniform_double(random_engine.rand_state)); ,
                        return random_util::one_to_zero(curand_uniform_double(random_engine.rand_state));
                )
        ))
#endif
        AMREX_IF_ON_HOST((
                amrex::ignore_unused(random_engine);
                return Random();
        ))
#endif
    }

    /**
    * \ingroup amrex_utilities
    * \brief Generate a psuedo-random real from a normal distribution
    *
    *  Generates one pseudorandom real number from a normal
    *  distribution with mean 'mean' and standard deviation 'stddev'.
    *
    */
    Real RandomNormal (Real mean, Real stddev);

    /// \ingroup amrex_utilities
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real RandomNormal (Real mean, Real stddev, RandomEngine const& random_engine)
    {
#if defined (__SYCL_DEVICE_ONLY__)
        mkl::rng::device::gaussian<Real> distr(mean, stddev);
        return mkl::rng::device::generate(distr, *random_engine.engine);
#else
#ifdef BL_USE_FLOAT
        AMREX_IF_ON_DEVICE((
                AMREX_HIP_OR_CUDA(
                        return stddev * hiprand_normal(random_engine.rand_state) + mean; ,
                        return stddev * curand_normal(random_engine.rand_state) + mean;
                )
        ))
#else
        AMREX_IF_ON_DEVICE((
                AMREX_HIP_OR_CUDA(
                        return stddev * hiprand_normal_double(random_engine.rand_state) + mean; ,
                        return stddev * curand_normal_double(random_engine.rand_state) + mean;
                )
        ))
#endif
        AMREX_IF_ON_HOST((
                amrex::ignore_unused(random_engine);
                return RandomNormal(mean, stddev);
        ))
#endif
    }

    /**
    * \ingroup amrex_utilities
    * \brief Generate a psuedo-random integer from a Poisson distribution
    *
    *  Generates one pseudorandom positive integer number
    *  extracted from a Poisson distribution, given the Real parameter lambda.
    *  The CPU version of this function relies on the standard Template Library
    *  The GPU version of this function relies on the cuRAND library
    *
    */
    unsigned int RandomPoisson (Real lambda);

    /// \ingroup amrex_utilities
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    unsigned int RandomPoisson (Real lambda, RandomEngine const& random_engine)
    {
#if defined (__SYCL_DEVICE_ONLY__)
        mkl::rng::device::poisson<unsigned int> distr(lambda);
        return mkl::rng::device::generate(distr, *random_engine.engine);
#else
        AMREX_IF_ON_DEVICE((
                AMREX_HIP_OR_CUDA(
                        return hiprand_poisson(random_engine.rand_state, lambda); ,
                        return curand_poisson(random_engine.rand_state, lambda);
                )
        ))
        AMREX_IF_ON_HOST((
                amrex::ignore_unused(random_engine);
                return RandomPoisson(lambda);
        ))
#endif
    }

    /// \cond DOXYGEN_IGNORE
    namespace random_util {

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    Real RandomGamma_alpha_ge_1 (Real alpha, Real beta, RandomEngine const& random_engine)
    {
        AMREX_ASSERT(alpha >= 1);
        AMREX_ASSERT(beta > 0);

        Real x, v, u;
        Real d = alpha - 1.0_rt / 3.0_rt;
        Real c = (1.0_rt / 3.0_rt) / std::sqrt(d);

        while (true) {
            do {
                x = amrex::RandomNormal(0.0_rt, 1.0_rt, random_engine);
                v = 1.0_rt + c * x;
            } while (v <= 0.0_rt);

            v = v * v * v;
            u = amrex::Random(random_engine);

            if (u < 1.0_rt - 0.0331_rt * x * x * x * x) {
                break;
            }

            if (std::log(u) < 0.5_rt * x * x + d * (1.0_rt - v + std::log(v))) {
                break;
            }
        }
        return beta * d * v;
    }
    }
    /// \endcond

    /**
    * \ingroup amrex_utilities
    * \brief Generate a psuedo-random floating point number from the Gamma distribution
    *
    *  Generates one real number (single or double)
    *  extracted from a Gamma distribution, given the Real parameters alpha and beta.
    *  alpha and beta must both be > 0.
    *  The CPU version of this function relies on the Standard Template Library
    *  The GPU version of this function relies is implemented in terms of Random
    *  and RandomNormal.
    *
    */
    Real RandomGamma (Real alpha, Real beta);

    /// \ingroup amrex_utilities
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real RandomGamma (Real alpha, Real beta, RandomEngine const& random_engine)
    {
        AMREX_ASSERT(alpha > 0);
        AMREX_ASSERT(beta > 0);

        AMREX_IF_ON_DEVICE((
        if (alpha < 1)
        {
            Real u = amrex::Random(random_engine);
            return amrex::random_util::RandomGamma_alpha_ge_1(1.0_rt + alpha, beta, random_engine) * std::pow(u, 1.0_rt / alpha);
        } else {
            return amrex::random_util::RandomGamma_alpha_ge_1(alpha, beta, random_engine);
        }
        ))

        AMREX_IF_ON_HOST((
            amrex::ignore_unused(random_engine);
            return RandomGamma(alpha, beta);
        ))
    }

    /**
    * \ingroup amrex_utilities
    * \brief Generates one pseudorandom unsigned integer which is
    *  uniformly distributed on [0,n-1]-interval for each call. For
    *  the special case of n = 0, it returns 0.
    *
    * The CPU version of this function uses C++11's mt19937.
    * The GPU version uses CURAND's XORWOW generator.
    */
    unsigned int Random_int (unsigned int n); // [0,n-1]

    /// \ingroup amrex_utilities
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    unsigned int Random_int (unsigned int n, RandomEngine const& random_engine)
    {
    if (n == 0) { return 0;}
#if defined(__SYCL_DEVICE_ONLY__)
        mkl::rng::device::uniform<unsigned int> distr(0,n);
        return mkl::rng::device::generate(distr, *random_engine.engine);
#else
        AMREX_IF_ON_DEVICE((
                unsigned int rand;
                constexpr unsigned int RAND_M = 4294967295; // 2**32-1
                do {
                AMREX_HIP_OR_CUDA( rand = hiprand(random_engine.rand_state);,
                                rand =  curand(random_engine.rand_state) );
                } while (rand >= (RAND_M - RAND_M % n));
                return rand % n;
        ))
        AMREX_IF_ON_HOST((
                amrex::ignore_unused(random_engine);
                return Random_int(n);
        ))
#endif
    }

    /**
    * \ingroup amrex_utilities
    * \brief Generates one pseudorandom unsigned long which is
    *  uniformly distributed on [0,n-1]-interval for each call. For
    *  the special case of n = 0, it returns 0.
    *
    * The CPU version of this function uses C++11's mt19937.
    * There is no GPU version.
    */
    ULong Random_long (ULong n); // [0,n-1]

    //! \ingroup amrex_utilities
    //! Fill random numbers from uniform distribution.  The interval is not
    //! uniform across backends the way amrex::Random is: it is [0,1) on the
    //! CPU, but the GPU paths hand back whatever the vendor bulk generator
    //! produces, nominally [0,1) for SYCL and (0,1] for CUDA and HIP.  Code
    //! that depends on either endpoint should guard the values it reads back,
    //! or draw with amrex::Random instead.

    //! \ingroup amrex_utilities
    //! Fill random numbers from normal distribution
    void FillRandomNormal (Real* p, Long N, Real mean, Real stddev);

    /// \cond DOXYGEN_IGNORE
    namespace detail {
        inline ULong DefaultGpuSeed () {
            return ParallelDescriptor::MyProc()*1234567ULL + 12345ULL;
        }
    }
    /// \endcond

    /** \brief Set the seed of the random number generator.
    * \ingroup amrex_utilities
    *
    *  There is also an entry point for Fortran callable as:
    *
    *  INTEGER seed
    *  call blutilinitrand(seed)
    *
    *  or
    *
    *  INTEGER seed
    *  call blinitrand(seed)
    */
    void InitRandom (ULong cpu_seed, int nprocs=ParallelDescriptor::NProcs(),
                     ULong gpu_seed = detail::DefaultGpuSeed());

    /// \ingroup amrex_utilities
    void ResetRandomSeed (ULong cpu_seed, ULong gpu_seed = detail::DefaultGpuSeed());

    /**
    * \ingroup amrex_utilities
    * \brief Save host random state.
    *
    * This function only serializes the host-side random number generators used
    * by the host random-number APIs. It does not save GPU random state used by
    * ParallelForRNG or GPU FillRandom/FillRandomNormal calls.
    */
    void SaveRandomState (std::ostream& os);

    /**
    * \ingroup amrex_utilities
    * \brief Restore host random state saved by SaveRandomState.
    *
    * This does not restore GPU random state used by ParallelForRNG or GPU
    * FillRandom/FillRandomNormal calls.
    */
    void RestoreRandomState (std::istream& is, int nthreads_old, int nstep_old);

    /**
    * \ingroup amrex_utilities
    * \brief Create a unique subset of random numbers from a pool
    *   of integers in the range [0, poolSize - 1]
    *   the set will be in the order they are found
    *   setSize must be <= poolSize
    *   uSet will be resized to setSize
    *   if you want all processors to have the same set,
    *   call this on one processor and broadcast the array
    */
    void UniqueRandomSubset (Vector<int> &uSet, int setSize, int poolSize,
                             bool printSet = false);

    void DeallocateRandomSeedDevArray ();
}

#endif
