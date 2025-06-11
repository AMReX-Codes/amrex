#ifndef BL_MEMPOOL_H
#define BL_MEMPOOL_H
#include <AMReX_Config.H>

#include <AMReX_GpuLaunch.H>
#include <AMReX_REAL.H>

#include <limits>
#include <type_traits>

extern "C" {
    void  amrex_mempool_init ();
    void  amrex_mempool_finalize ();
    void* amrex_mempool_alloc (size_t n);
    void  amrex_mempool_free (void* p);
    void  amrex_mempool_get_stats (int& mp_min, int& mp_max, int& mp_tot);  //!< min, max & tot in MB
    void  amrex_real_array_init (amrex_real* p, size_t nelems);
    void  amrex_array_init_snan (amrex_real* p, size_t nelems);
}

namespace amrex {
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
