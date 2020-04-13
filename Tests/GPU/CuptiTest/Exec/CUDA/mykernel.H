#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>

using clock_value_t = long long;

AMREX_GPU_DEVICE
void deviceSleep (clock_value_t sleep_cycles)
{
    clock_value_t start = clock64();
    clock_value_t cycles_elapsed;
    do { cycles_elapsed = clock64() - start; } 
    while (cycles_elapsed < sleep_cycles);
}

AMREX_GPU_DEVICE
void deviceNestedSleep (clock_value_t sleep_cycles, int level)
{
    if (level==0) {
        return;
    }
    clock_value_t start = clock64();
    clock_value_t cycles_elapsed;
    do {cycles_elapsed = clock64() - start;} 
    while (cycles_elapsed < sleep_cycles);

    deviceNestedSleep(sleep_cycles, level - 1);
}

AMREX_GPU_GLOBAL
void deviceSleep ()
{
    // Sleep for some number of cycles
    int cycles;
    cycles = (int) (1e8);
    deviceSleep(cycles);
}

#endif
