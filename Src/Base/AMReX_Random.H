#ifndef AMREX_RAND_H
#define AMREX_RAND_H

#include <AMReX_GpuQualifiers.H>
#include <AMReX_ParallelDescriptor.H>

// xxxxx DPCPP todo

namespace amrex
{
    /**
    * \brief Generate a psuedo-random double using C++11's mt19937.
    *
    *  Generates one pseudorandom real number (double) from a normal
    *  distribution with mean 'mean' and standard deviation 'stddev'.
    *
    */
    AMREX_GPU_EXTERNAL AMREX_GPU_HOST_DEVICE Real RandomNormal (Real mean, Real stddev);

    /**
    * \brief Generate a psuedo-random double using C++11's mt19937.
    *
    *  Generates one pseudorandom real number (double) from a uniform
    *  distribution between 0.0 and 1.0 (0.0 included, 1.0 excluded)
    *
    */
    AMREX_GPU_EXTERNAL AMREX_GPU_HOST_DEVICE Real Random ();

    /**
    * \brief Generate a psuedo-random integer from a Poisson distribution
    *
    *  Generates one pseudorandom positive integer number (double)
    *  extracted from a Poisson distribution, given the Real parameter lambda.
    *  The CPU version of this function relies on the standard Template Library
    *  The GPU version of this function relies on the cuRAND library
    *
    */
    AMREX_GPU_EXTERNAL AMREX_GPU_HOST_DEVICE unsigned int RandomPoisson (Real lambda);



#ifdef AMREX_USE_GPU
    // Locking mechanism functions for locking and unlocking
    AMREX_GPU_DEVICE int get_state (int tid);
    AMREX_GPU_DEVICE void free_state (int tid);
#endif


    /**
    * \brief Generates one pseudorandom unsigned integer which is
    *  uniformly distributed on [0,n-1]-interval for each call.
    *
    * The CPU version of this function uses C++11's mt19937.
    * The GPU version uses CURAND's XORWOW generator.
    */
    AMREX_GPU_EXTERNAL AMREX_GPU_HOST_DEVICE unsigned int Random_int (unsigned int n); // [0,n-1]

    /**
    * \brief Generates one pseudorandom unsigned long which is
    *  uniformly distributed on [0,n-1]-interval for each call.
    *
    * The CPU version of this function uses C++11's mt19937.
    * There is no GPU version.
    */
    AMREX_GPU_HOST ULong Random_long (ULong n); // [0,n-1]

    /** \brief Set the seed of the random number generator.
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
    void InitRandom (ULong seed, int nprocs=ParallelDescriptor::NProcs());

    /**
    * \brief Resize seed array and copy address to global symbol
    */
    void ResizeRandomSeed (int N);
    void ResetRandomSeed (ULong seed);

    /**
    * \brief Set the seeds of the random number generator for each CUDA thread .
    *  An in-built CUDA library, curand_init is used here.
    *  The pseusorandom sequence currently implemented is obtained
    *  from the XORWOW and MRG32k3a generator
    */
    void InitRandSeedOnDevice (int N);

    /**
    * \brief Save and restore random state.
    *
    */
    void SaveRandomState (std::ostream& os);

    void RestoreRandomState (std::istream& is, int nthreads_old, int nstep_old);

    /**
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

    /**
    * \brief Fill binCounts with the number of items for each bin.
    *   The remainder items are randomly assigned.
    */
    void NItemsPerBin (int totalItems, Vector<int> &binCounts);

    /**
    * \brief Set the seeds of the random number generator for each CUDA thread .
    *  An in-built CUDA library, curand_init is used here.
    *  The pseusorandom sequence currently implemented is obtained
    *  from the XORWOW and MRG32k3a generator
    */
    void InitRandSeedOnDevice (int N);

    void DeallocateRandomSeedDevArray ();
}

#endif
