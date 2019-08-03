#ifndef AMREX_PARTITION_H_
#define AMREX_PARTITION_H_

#include <AMReX_Gpu.H>
#include <AMReX_Scan.H>

namespace amrex {

#ifdef AMREX_USE_GPU
template <typename T, typename F>
int Partition (Gpu::DeviceVector<T> x, F && f)
{
    return 0;
}
#endif

}

#endif
