#ifndef AMREX_CUDA_ASYNC_FAB_H_
#define AMREX_CUDA_ASYNC_FAB_H_

#include <memory>
#include <AMReX_Box.H>

namespace amrex {

class FArrayBox;

namespace Cuda {

class AsyncFabImpl;

class AsyncFab
{
public:

    AsyncFab ();
    AsyncFab (Box const& bx, int ncomp);
    AsyncFab (FArrayBox& a_fab);
    AsyncFab (FArrayBox& a_fab, Box const& bx, int comp);
#ifdef AMREX_USE_GPU
    AMREX_GPU_DEVICE
    AsyncFab (AsyncFab const& rhs) : m_impl(nullptr), m_fab(rhs.m_fab) {}
#else
    AsyncFab (AsyncFab const& rhs) = delete;
#endif

    AsyncFab (AsyncFab&& rhs) : m_impl(rhs.m_impl), m_fab(rhs.m_fab) {
        rhs.m_impl = nullptr;  rhs.m_fab = nullptr;
    }

    AMREX_GPU_HOST_DEVICE
    ~AsyncFab () {
#ifndef __CUDA_ARCH__
        clear();
#endif
    }

    void operator= (AsyncFab const&) = delete;

// Not allowed because it contains cuda api calls   void resize (Box const& bx, int ncomp);
    void clear ();
    AMREX_GPU_HOST_DEVICE FArrayBox* fabPtr () const noexcept { return m_fab; }
    AMREX_GPU_HOST_DEVICE FArrayBox& fab () const noexcept { return *m_fab; }
    FArrayBox& hostFab () const noexcept;

    Array4<Real const> array () const noexcept;
    Array4<Real> array () noexcept;

    static void Initialize ();
    static void Finalize ();

private:
    AsyncFabImpl* m_impl = nullptr;
    FArrayBox* m_fab = nullptr;
};

}

using Cuda::AsyncFab;

}

#endif
