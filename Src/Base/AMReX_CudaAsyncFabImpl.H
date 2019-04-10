#ifndef AMREX_CUDA_ASYNC_FAB_IMPL_H_
#define AMREX_CUDA_ASYNC_FAB_IMPL_H_

#include <memory>
#include <AMReX_FArrayBox.H>

namespace amrex {
namespace Cuda {

class AsyncFabImpl
{
public:

    AsyncFabImpl ();
    AsyncFabImpl (Box const& bx, int ncomp);
    AsyncFabImpl (FArrayBox& a_fab);
    AsyncFabImpl (FArrayBox& a_fab, Box const& bx, int ncomp);
    ~AsyncFabImpl ();
    AsyncFabImpl (AsyncFabImpl const&) = delete;
    void operator= (AsyncFabImpl const&) = delete;

    FArrayBox* fabPtr () noexcept;
    FArrayBox& hostFab () noexcept;

    static void Initialize ();
    static void Finalize ();

private:

#ifdef AMREX_USE_CUDA

#ifdef AMREX_FAB_IS_PINNED

    std::unique_ptr<FArrayBox> m_gpu_fab;

#else

    FArrayBox                  m_cpu_fab;
    std::unique_ptr<FArrayBox> m_gpu_fab;
    BaseFabData<Real>          m_cpu_fab_data;

    void copy_htod ();

#endif

#else
    // cpu version
    FArrayBox                  m_cpu_fab;
    FArrayBox*                 m_cpu_fab_alias = nullptr;
#endif

};

}
}

#endif
