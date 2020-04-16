#include <cstring>
#include <cstdlib>

#include <AMReX_BaseFab.H>
#include <AMReX_BLFort.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

namespace amrex {

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
    return r;
#else
    return private_total_bytes_allocated_in_fabs;
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
    return r;
#else
    return private_total_bytes_allocated_in_fabs_hwm;
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
    return r;
#else
    return private_total_cells_allocated_in_fabs;
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
    return r;
#else
    return private_total_cells_allocated_in_fabs_hwm;
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
}

void
update_fab_stats (Long n, Long s, size_t szt) noexcept
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
}

}
