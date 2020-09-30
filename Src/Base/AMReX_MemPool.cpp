#include <iostream>
#include <limits>
#include <algorithm>
#include <new>
#include <memory>
#include <cstring>
#include <cstdint>

#include <AMReX_CArena.H>
#include <AMReX_MemPool.H>
#include <AMReX_Vector.H>
#include <AMReX_OpenMP.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <AMReX_ParmParse.H>

using namespace amrex;

namespace
{
    static Vector<std::unique_ptr<CArena> > the_memory_pool;
#if defined(AMREX_TESTING) || defined(AMREX_DEBUG)
    static int init_snan = 1;
#else
    static int init_snan = 0;
#endif
    static bool initialized = false;
}

extern "C" {

void amrex_mempool_init ()
{
    if (!initialized)
    {
	initialized = true;

        ParmParse pp("fab");
	pp.query("init_snan", init_snan);

	int nthreads = OpenMP::get_max_threads();

	the_memory_pool.resize(nthreads);
	for (int i=0; i<nthreads; ++i) {
// xxxxx HIP FIX THIS - Default Arena w/o managed?
// Default arena is currently Device on HIP where there is no managed option.
// Need to adjust to CPU specifically in that case.
#ifdef AMREX_USE_HIP
            the_memory_pool[i].reset(new CArena(0, ArenaInfo().SetCpuMemory()));
#else
            the_memory_pool[i].reset(new CArena);
#endif
	}

#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
	{
	    size_t N = 1024*1024*sizeof(double);
	    void *p = amrex_mempool_alloc(N);
	    memset(p, 0, N);
	    amrex_mempool_free(p);
	}

#ifdef AMREX_MEM_PROFILING
	MemProfiler::add("MemPool", std::function<MemProfiler::MemInfo()>
			 ([] () -> MemProfiler::MemInfo {
			     int MB_min, MB_max, MB_tot;
			     amrex_mempool_get_stats(MB_min, MB_max, MB_tot);
			     Long b = MB_tot * (1024L*1024L);
			     return {b, b};
			 }));
#endif
    }
}

void amrex_mempool_finalize ()
{
    initialized = false;
    the_memory_pool.clear();
}

void* amrex_mempool_alloc (size_t nbytes)
{
  int tid = OpenMP::get_thread_num();
  return the_memory_pool[tid]->alloc(nbytes);
}

void amrex_mempool_free (void* p) 
{
  int tid = OpenMP::get_thread_num();
  the_memory_pool[tid]->free(p);
}

void amrex_mempool_get_stats (int& mp_min, int& mp_max, int& mp_tot) // min, max & tot in MB
{
  size_t hsu_min=std::numeric_limits<size_t>::max();
  size_t hsu_max=0;
  size_t hsu_tot=0;
  for (const auto& mp : the_memory_pool) {
    size_t hsu = mp->heap_space_used();
    hsu_min = std::min(hsu, hsu_min);
    hsu_max = std::max(hsu, hsu_max);
    hsu_tot += hsu;
  }
  mp_min = hsu_min/(1024*1024);
  mp_max = hsu_max/(1024*1024);
  mp_tot = hsu_tot/(1024*1024);
}

void amrex_real_array_init (Real* p, size_t nelems)
{
    if (init_snan) amrex_array_init_snan(p, nelems);
}

void amrex_array_init_snan (Real* p, size_t nelems)
{
#ifdef BL_USE_DOUBLE

#ifdef UINT64_MAX
    const uint64_t snan = UINT64_C(0x7ff0000080000001);
    static_assert(sizeof(double) == sizeof(uint64_t), "MemPool: sizeof double != sizeof uint64_t");
    for (size_t i = 0; i < nelems; ++i) {
        std::memcpy(p++, &snan, sizeof(double));
    }
#endif

#else

#ifdef UINT32_MAX
    const uint32_t snan = UINT32_C(0x7fa00000);
    static_assert(sizeof(float) == sizeof(uint32_t), "MemPool: sizeof float != sizeof uint32_t");
    for (size_t i = 0; i < nelems; ++i) {
        std::memcpy(p++, &snan, sizeof(float));
    }
#endif

#endif
}
}
