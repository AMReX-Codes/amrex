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
    ~AsyncFab ();
    AsyncFab (AsyncFab const&) = delete;
    void operator= (AsyncFab const&) = delete;

// Not allowed because it contains cuda api calls   void resize (Box const& bx, int ncomp);
    void clear ();
    FArrayBox* fabPtr ();

    static void Initialize ();
    static void Finalize ();

private:
    std::unique_ptr<AsyncFabImpl> m_impl;
};

}
}

#endif
