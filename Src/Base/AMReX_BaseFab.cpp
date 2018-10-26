#include <cstring>
#include <cstdlib>

#include <AMReX_BaseFab.H>
#include <AMReX_BArena.H>
#include <AMReX_CArena.H>

#include <AMReX_BLFort.H>

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

namespace amrex {

long private_total_bytes_allocated_in_fabs     = 0L;
long private_total_bytes_allocated_in_fabs_hwm = 0L;
long private_total_cells_allocated_in_fabs     = 0L;
long private_total_cells_allocated_in_fabs_hwm = 0L;

int BF_init::m_cnt = 0;

namespace
{
    static bool basefab_initialized = false;

    Arena* the_arena = nullptr;
    Arena* the_cuda_arena = nullptr;
    Arena* the_managed_arena = nullptr;
    Arena* the_pinned_arena = nullptr;
}

BF_init::BF_init ()
{
    if (m_cnt++ == 0)
    {
        BL_ASSERT(the_arena == nullptr);
        BL_ASSERT(the_cuda_arena == nullptr);
        BL_ASSERT(the_managed_arena == nullptr);
        BL_ASSERT(the_pinned_arena == nullptr);

#if defined(BL_COALESCE_FABS)
        the_arena = new CArena;
#else
        the_arena = new BArena;
#endif

#ifdef AMREX_USE_CUDA
        the_arena->SetPreferred();
#endif

#if AMREX_USE_CUDA
        the_cuda_arena = new CArena;
        the_cuda_arena->SetDeviceMemory();
#else
        the_cuda_arena = new BArena;
#endif

#if defined(AMREX_USE_CUDA)
        the_managed_arena = new CArena;
#else
        the_managed_arena = new BArena;
#endif

#if defined(AMREX_USE_CUDA)
        const std::size_t hunk_size = 64 * 1024;
        the_pinned_arena = new CArena(hunk_size);
        the_pinned_arena->SetHostAlloc();
#else
        the_pinned_arena = new BArena;
#endif
    }
}

BF_init::~BF_init ()
{
    if (--m_cnt == 0) {
        delete the_arena;
        the_arena = nullptr;

        delete the_cuda_arena;
        the_cuda_arena = nullptr;

        delete the_managed_arena;
        the_managed_arena = nullptr;

        delete the_pinned_arena;
        the_pinned_arena = nullptr;
    }
}

Arena*
The_Arena ()
{
    BL_ASSERT(the_arena != nullptr);
    return the_arena;
}

Arena*
The_Cuda_Arena ()
{
    BL_ASSERT(the_cuda_arena != nullptr);
    return the_cuda_arena;
}

Arena*
The_Managed_Arena ()
{
    BL_ASSERT(the_managed_arena != nullptr);
    return the_managed_arena;
}

Arena*
The_Pinned_Arena ()
{
    BL_ASSERT(the_pinned_arena != nullptr);
    return the_pinned_arena;
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
}


long 
TotalBytesAllocatedInFabs()
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
TotalBytesAllocatedInFabsHWM()
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
TotalCellsAllocatedInFabs()
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
TotalCellsAllocatedInFabsHWM()
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
ResetTotalBytesAllocatedInFabsHWM()
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	private_total_bytes_allocated_in_fabs_hwm = 0;
    }
}

void
update_fab_stats (long n, long s, size_t szt)
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
