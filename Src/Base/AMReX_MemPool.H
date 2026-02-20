#ifndef BL_MEMPOOL_H
#define BL_MEMPOOL_H
#include <AMReX_Config.H>

/**
 * \file AMReX_MemPool.H
 * \brief C-linkage memory pool interface and the \c fill_snan utility.
 */

#include <AMReX_GpuLaunch.H>
#include <AMReX_REAL.H>

#include <limits>
#include <type_traits>

extern "C" {
    //! Initialize the AMReX memory pool (called at AMReX initialization).
    void  amrex_mempool_init ();
    //! Finalize the AMReX memory pool (called at AMReX finalization).
    void  amrex_mempool_finalize ();
    //! Allocate \p n bytes from the memory pool.
    void* amrex_mempool_alloc (size_t n);
    //! Free memory previously allocated with \c amrex_mempool_alloc.
    void  amrex_mempool_free (void* p);
    void  amrex_mempool_get_stats (int& mp_min, int& mp_max, int& mp_tot);  //!< min, max & tot in MB
    //! Initialize \p nelems Real values at \p p to zero.
    void  amrex_real_array_init (amrex_real* p, size_t nelems);
    //! Fill \p nelems Real values at \p p with signaling NaN.
    void  amrex_array_init_snan (amrex_real* p, size_t nelems);
}

namespace amrex {
    /**
     * \ingroup amrex_memory
     * \brief Fill \p nelems floating-point values at \p p with signaling NaN.
     * \tparam run_on Execute on \c RunOn::Device (GPU) or \c RunOn::Host (CPU).
     */
    template <RunOn run_on, typename T,
              std::enable_if_t<std::is_same_v<T,double> || std::is_same_v<T,float>, int> FOO = 0>
    void fill_snan (T* p, std::size_t nelems)
    {
        if (p == nullptr || nelems == 0) { return; }
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion() && run_on == RunOn::Device) {
            amrex::ParallelFor(nelems, [=] AMREX_GPU_DEVICE (Long i) noexcept
            {
                p[i] = std::numeric_limits<T>::signaling_NaN();
            });
        } else
#endif
        {
            if constexpr (std::is_same_v<T,float>) {
#ifdef UINT32_MAX
                const uint32_t snan = UINT32_C(0x7fa00000);
                static_assert(sizeof(float) == sizeof(uint32_t), "MemPool: sizeof float != sizeof uint32_t");
                for (size_t i = 0; i < nelems; ++i) {
                    std::memcpy(p++, &snan, sizeof(float));
                }
#endif
            } else if constexpr (std::is_same_v<T,double>) {
#ifdef UINT64_MAX
                const uint64_t snan = UINT64_C(0x7ff0000080000001);
                static_assert(sizeof(double) == sizeof(uint64_t), "MemPool: sizeof double != sizeof uint64_t");
                for (size_t i = 0; i < nelems; ++i) {
                    std::memcpy(p++, &snan, sizeof(double));
                }
#endif
            }
        }
    }
}

#endif
