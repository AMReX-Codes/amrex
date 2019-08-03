#ifndef AMREX_PARTITION_H_
#define AMREX_PARTITION_H_

#include <AMReX_Gpu.H>

namespace amrex {

#ifdef AMREX_USE_GPU
template <typename T, typename F>
int Parition (Gpu::DeviceDevice<T> x, F && f)
{
    return 0;
}
#endif

}

#endif
