#include <cstring>
#include <cstdlib>

#include <AMReX_BaseFab.H>
#include <AMReX_BLFort.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

namespace amrex {

std::atomic<Long> atomic_total_bytes_allocated_in_fabs     {0L};
std::atomic<Long> atomic_total_bytes_allocated_in_fabs_hwm {0L};
std::atomic<Long> atomic_total_cells_allocated_in_fabs     {0L};
std::atomic<Long> atomic_total_cells_allocated_in_fabs_hwm {0L};
Long private_total_bytes_allocated_in_fabs     = 0L;
Long private_total_bytes_allocated_in_fabs_hwm = 0L;
Long private_total_cells_allocated_in_fabs     = 0L;
Long private_total_cells_allocated_in_fabs_hwm = 0L;

namespace
{
    static bool basefab_initialized = false;
}

void
BaseFab_Initialize ()
{
    if (!basefab_initialized)
    {
        basefab_initialized = true;

#ifdef _OPENMP
#pragma omp parallel
        {
            amrex::private_total_bytes_allocated_in_fabs     = 0;
            amrex::private_total_bytes_allocated_in_fabs_hwm = 0;
            amrex::private_total_cells_allocated_in_fabs     = 0;
            amrex::private_total_cells_allocated_in_fabs_hwm = 0;
        }
#endif

#ifdef AMREX_MEM_PROFILING
        MemProfiler::add("Fab", std::function<MemProfiler::MemInfo()>
                         ([] () -> MemProfiler::MemInfo {
                             return {amrex::TotalBytesAllocatedInFabs(),
                                     amrex::TotalBytesAllocatedInFabsHWM()};
                         }));
#endif
    }

    amrex::ExecOnFinalize(amrex::BaseFab_Finalize);
}

void
BaseFab_Finalize()
{
    basefab_initialized = false;
}


Long
TotalBytesAllocatedInFabs () noexcept
{
#ifdef _OPENMP
    Long r=0;
#pragma omp parallel reduction(+:r)
    {
        r += private_total_bytes_allocated_in_fabs;
    }
    return r
        + atomic_total_bytes_allocated_in_fabs.load(std::memory_order_relaxed);
#else
    return private_total_bytes_allocated_in_fabs
        + atomic_total_bytes_allocated_in_fabs.load(std::memory_order_relaxed);
#endif
}

Long
TotalBytesAllocatedInFabsHWM () noexcept
{
#ifdef _OPENMP
    Long r=0;
#pragma omp parallel reduction(+:r)
    {
        r += private_total_bytes_allocated_in_fabs_hwm;
    }
    return r
        + atomic_total_bytes_allocated_in_fabs_hwm.load(std::memory_order_relaxed);
#else
    return private_total_bytes_allocated_in_fabs_hwm
        + atomic_total_bytes_allocated_in_fabs_hwm.load(std::memory_order_relaxed);
#endif
}

Long
TotalCellsAllocatedInFabs () noexcept
{
#ifdef _OPENMP
    Long r=0;
#pragma omp parallel reduction(+:r)
    {
        r += private_total_cells_allocated_in_fabs;
    }
    return r
        + atomic_total_cells_allocated_in_fabs.load(std::memory_order_relaxed);
#else
    return private_total_cells_allocated_in_fabs
        + atomic_total_cells_allocated_in_fabs.load(std::memory_order_relaxed);
#endif
}

Long
TotalCellsAllocatedInFabsHWM () noexcept
{
#ifdef _OPENMP
    Long r=0;
#pragma omp parallel reduction(+:r)
    {
        r += private_total_cells_allocated_in_fabs_hwm;
    }
    return r
        + atomic_total_cells_allocated_in_fabs_hwm.load(std::memory_order_relaxed);
#else
    return private_total_cells_allocated_in_fabs_hwm
        + atomic_total_cells_allocated_in_fabs_hwm.load(std::memory_order_relaxed);
#endif
}

void 
ResetTotalBytesAllocatedInFabsHWM () noexcept
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        private_total_bytes_allocated_in_fabs_hwm = 0;
    }
    atomic_total_bytes_allocated_in_fabs_hwm.store(0,std::memory_order_relaxed);
}

void
update_fab_stats (Long n, Long s, size_t szt) noexcept
{
#ifdef _OPENMP
    if (omp_in_parallel())
    {
        Long tst = s*szt;
        amrex::private_total_bytes_allocated_in_fabs += tst;
        amrex::private_total_bytes_allocated_in_fabs_hwm
            = std::max(amrex::private_total_bytes_allocated_in_fabs_hwm,
                       amrex::private_total_bytes_allocated_in_fabs);

        if(szt == sizeof(Real)) {
            amrex::private_total_cells_allocated_in_fabs += n;
            amrex::private_total_cells_allocated_in_fabs_hwm
                = std::max(amrex::private_total_cells_allocated_in_fabs_hwm,
                           amrex::private_total_cells_allocated_in_fabs);
        }
    } else
#endif
    {
        Long tst = s*szt;
        Long old_bytes = amrex::atomic_total_bytes_allocated_in_fabs.fetch_add
            (tst,std::memory_order_relaxed);
        Long new_bytes = old_bytes + tst;
        Long prev_bytes_hwm = amrex::atomic_total_bytes_allocated_in_fabs_hwm.load
            (std::memory_order_relaxed);
        while (prev_bytes_hwm < new_bytes) {
            if (amrex::atomic_total_bytes_allocated_in_fabs_hwm.compare_exchange_weak
                (prev_bytes_hwm, new_bytes, std::memory_order_release, std::memory_order_relaxed)) {
                break;
            }
        }

        if(szt == sizeof(Real)) {
            Long old_cells = amrex::atomic_total_cells_allocated_in_fabs.fetch_add
                (n,std::memory_order_relaxed);
            Long new_cells = old_cells + n;
            Long prev_cells_hwm = amrex::atomic_total_cells_allocated_in_fabs_hwm.load
                (std::memory_order_relaxed);
            while (prev_cells_hwm < new_cells) {
                if (amrex::atomic_total_cells_allocated_in_fabs_hwm.compare_exchange_weak
                    (prev_cells_hwm, new_cells, std::memory_order_release, std::memory_order_relaxed)) {
                    break;
                }
            }
        }
    }
}

}
