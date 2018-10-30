#ifndef AMREX_CUDA_FAB_H_
#define AMREX_CUDA_FAB_H_

#include <AMReX_Box.H>

namespace amrex {

class FArrayBox;

namespace Cuda {

class DeviceFab
{
public:

    DeviceFab ();
    DeviceFab (Box const& bx, int ncomp);
    ~DeviceFab ();
    DeviceFab (DeviceFab const&) = delete;
    void operator= (DeviceFab const&) = delete;

private:
    FArrayBox* m_fab = nullptr;
};

}
}

#endif
