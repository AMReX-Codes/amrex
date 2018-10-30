#ifndef AMREX_CUDA_FAB_H_
#define AMREX_CUDA_FAB_H_

#include <memory>
#include <AMReX_Box.H>

namespace amrex {

class FArrayBox;

namespace Cuda {

class DeviceFabImpl;

class DeviceFab
{
public:

    DeviceFab ();
    DeviceFab (Box const& bx, int ncomp);
    ~DeviceFab ();
    DeviceFab (DeviceFab const&) = delete;
    void operator= (DeviceFab const&) = delete;

    void resize (Box const& bx, int ncomp);
    FArrayBox* fabPtr ();

    static void Initialize ();
    static void Finalize ();

private:
    std::unique_ptr<DeviceFabImpl> m_impl;
};

}
}

#endif
