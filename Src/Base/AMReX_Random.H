#ifndef AMREX_RAND_H
#define AMREX_RAND_H

#include <AMReX_GpuQualifiers.H>
#include <AMReX_ParallelDescriptor.H>

namespace amrex
{
    /**
    * \brief Generate a psuedo-random double using C++11's mt19937.
    *
    *  Generates one pseudorandom real number (double) from a normal
    *  distribution with mean 'mean' and standard deviation 'stddev'.
    *
    */
    AMREX_GPU_HOST_DEVICE Real RandomNormal (Real mean, Real stddev);

    /**
    * \brief Generate a psuedo-random double using C++11's mt19937.
    *
    *  Generates one pseudorandom real number (double) which is
    *  uniformly distributed on [0,1)-interval for each call.
    *
    */
    AMREX_GPU_HOST_DEVICE Real Random (); // [0,1)



#ifdef AMREX_USE_CUDA
    // Locking mechanism functions for locking and unlocking
    AMREX_GPU_DEVICE int get_state(int tid);
    AMREX_GPU_DEVICE void free_state(int tid);
#endif


    /**
    * \brief Generates one pseudorandom unsigned long integer which is
    *  uniformly distributed on [0,n-1]-interval for each call.
    *
    * The CPU version of this function uses C++11's mt19937.
    * The GPU version uses CURAND's XORWOW generator.
    */
#ifdef __CUDA_ARCH__    
    AMREX_GPU_DEVICE
    unsigned int Random_int(unsigned int n); // [0,n-1]
#else
    AMREX_GPU_HOST
    unsigned long Random_int(unsigned long n); // [0,n-1]    
#endif

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
    void InitRandom (unsigned long seed, int nprocs=ParallelDescriptor::NProcs());

    /**
    * \brief Resize seed array and copy address to global symbol
    */
    void ResizeRandomSeed(int N);
    void ResetRandomSeed(unsigned long seed);

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
    void SaveRandomState(std::ostream& os);

    void RestoreRandomState(std::istream& is, int nthreads_old, int nstep_old);

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
