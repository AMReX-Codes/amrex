#ifdef _OPENMP
#include <omp.h>
#endif

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

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <AMReX_ParmParse.H>

#ifdef USE_PERILLA_PTHREADS
#include <WorkerThread.H>
#endif

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

	int nthreads = 1;

#ifdef _OPENMP
	nthreads = omp_get_max_threads();
#endif


#ifdef USE_PERILLA_PTHREADS
#ifdef _OPENMP
	//Just in case Perilla thread spawns multiple OMP threads
        nthreads *= perilla::nThreads();
#else
	nthreads = perilla::nThreads();
#endif
#endif

	the_memory_pool.resize(nthreads);
	for (int i=0; i<nthreads; ++i) {
	    the_memory_pool[i].reset(new CArena);
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

#ifdef BL_MEM_PROFILING
	MemProfiler::add("MemPool", std::function<MemProfiler::MemInfo()>
			 ([] () -> MemProfiler::MemInfo {
			     int MB_min, MB_max, MB_tot;
			     amrex_mempool_get_stats(MB_min, MB_max, MB_tot);
			     long b = MB_tot * (1024L*1024L);
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
  int tid=0;

#ifdef _OPENMP
  tid = omp_get_thread_num();
#endif

#ifdef USE_PERILLA_PTHREADS
#ifdef _OPENMP
  tid = perilla::tid()*omp_get_max_threads()+tid;
#else
  tid = perilla::tid();
#endif
#endif
  return the_memory_pool[tid]->alloc(nbytes);
}

void amrex_mempool_free (void* p) 
{
  int tid=0;

#ifdef _OPENMP
  tid = omp_get_thread_num();
#endif

#ifdef USE_PERILLA_PTHREADS
#ifdef _OPENMP
  tid = perilla::tid()*omp_get_max_threads()+tid;
#else
  tid = perilla::tid();
#endif
#endif

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
#ifdef UINT64_MAX
    const uint64_t snan = UINT64_C(0x7ff0000080000001);
#else
    static_assert(sizeof(double) == sizeof(long long), "MemPool: sizeof double != sizeof long long");
    const long long snan = 0x7ff0000080000001LL;
#endif
    for (size_t i = 0; i < nelems; ++i) {
        std::memcpy(p++, &snan, sizeof(double));
    }
}
}
