#include <cstring>
#include <cstdlib>

#include <AMReX_BaseFab.H>
#include <AMReX_BLFort.H>

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

namespace amrex {

long private_total_bytes_allocated_in_fabs     = 0L;
long private_total_bytes_allocated_in_fabs_hwm = 0L;
long private_total_cells_allocated_in_fabs     = 0L;
long private_total_cells_allocated_in_fabs_hwm = 0L;

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

#ifdef AMREX_USE_GPU
        makeFabPoolAllocator();
#endif

#ifdef _OPENMP
#pragma omp parallel
        {
            amrex::private_total_bytes_allocated_in_fabs     = 0;
            amrex::private_total_bytes_allocated_in_fabs_hwm = 0;
            amrex::private_total_cells_allocated_in_fabs     = 0;
            amrex::private_total_cells_allocated_in_fabs_hwm = 0;
        }
#endif

#ifdef BL_MEM_PROFILING
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
#ifdef AMREX_USE_GPU
    destroyFabPoolAllocator();
#endif
}


long 
TotalBytesAllocatedInFabs () noexcept
{
#ifdef _OPENMP
    long r=0;
#pragma omp parallel reduction(+:r)
    {
	r += private_total_bytes_allocated_in_fabs;
    }
    return r;
#else
    return private_total_bytes_allocated_in_fabs;
#endif
}

long 
TotalBytesAllocatedInFabsHWM () noexcept
{
#ifdef _OPENMP
    long r=0;
#pragma omp parallel reduction(+:r)
    {
	r += private_total_bytes_allocated_in_fabs_hwm;
    }
    return r;
#else
    return private_total_bytes_allocated_in_fabs_hwm;
#endif
}

long 
TotalCellsAllocatedInFabs () noexcept
{
#ifdef _OPENMP
    long r=0;
#pragma omp parallel reduction(+:r)
    {
	r += private_total_cells_allocated_in_fabs;
    }
    return r;
#else
    return private_total_cells_allocated_in_fabs;
#endif
}

long 
TotalCellsAllocatedInFabsHWM () noexcept
{
#ifdef _OPENMP
    long r=0;
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
update_fab_stats (long n, long s, size_t szt) noexcept
{
    long tst = s*szt;
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
