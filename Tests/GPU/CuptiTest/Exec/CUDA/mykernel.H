#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void init_mf (int i, int j, int k,
               amrex::Array4<amrex::Real> const& mf)
{
  mf(i,j,k) = 0.;
}

using clock_value_t = long long;

__device__ void d_sleep(clock_value_t sleep_cycles)
{
  clock_value_t start = clock64();
  clock_value_t cycles_elapsed;
  do { cycles_elapsed = clock64() - start; } 
  while (cycles_elapsed < sleep_cycles);
}

__global__
void randomSleep ()
{
  // Sleep for random number of cycles
  int cycles;
  cycles = (int) (100000*amrex::Random());
  d_sleep(cycles);
}

#endif
