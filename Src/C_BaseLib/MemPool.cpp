#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <limits>
#include <algorithm>
#include <new>
#include <cstring>

#include <CArena.H>
#include <PArray.H>
#include <MemPool.H>

#ifdef BL_MEM_PROFILING
#include <MemProfiler.H>
#endif

#ifndef FORTRAN_BOXLIB
#include <ParmParse.H>
#endif

namespace
{
    static PArray<CArena> the_memory_pool;
#if defined(BL_TESTING) || defined(DEBUG)
    static int init_snan = 1;
#else
    static int init_snan = 0;
#endif
}

extern "C" {

void mempool_init()
{
    static bool initialized = false;
    if (!initialized)
    {
	initialized = true;

#ifndef FORTRAN_BOXLIB
        ParmParse pp("fab");
	pp.query("init_snan", init_snan);
#endif

#ifdef _OPENMP
	int nthreads = omp_get_max_threads();
#else
	int nthreads = 1;
#endif
	the_memory_pool.resize(nthreads, PArrayManage);
	for (int i=0; i<nthreads; ++i) {
	    the_memory_pool.set(i, new CArena());
	}
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
	    size_t N = 1024*1024*sizeof(double);
	    void *p = mempool_alloc(N);
	    memset(p, 0, N);
	    mempool_free(p);
	}

#ifdef BL_MEM_PROFILING
	MemProfiler::add("MemPool", std::function<MemProfiler::MemInfo()>
			 ([] () -> MemProfiler::MemInfo {
			     int MB_min, MB_max, MB_tot;
			     mempool_get_stats(MB_min, MB_max, MB_tot);
			     long b = MB_tot * (1024L*1024L);
			     return {b, b};
			 }));
#endif
    }
}

void* mempool_alloc (size_t nbytes)
{
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  return the_memory_pool[tid].alloc(nbytes);
}

void mempool_free (void* p) 
{
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  the_memory_pool[tid].free(p);
}

void mempool_get_stats (int& mp_min, int& mp_max, int& mp_tot) // min, max & tot in MB
{
  size_t hsu_min=std::numeric_limits<size_t>::max();
  size_t hsu_max=0;
  size_t hsu_tot=0;
  for (int i=0; i<the_memory_pool.size(); ++i) {
    size_t hsu = the_memory_pool[i].heap_space_used();
    hsu_min = std::min(hsu, hsu_min);
    hsu_max = std::max(hsu, hsu_max);
    hsu_tot += hsu;
  }
  mp_min = hsu_min/(1024*1024);
  mp_max = hsu_max/(1024*1024);
  mp_tot = hsu_tot/(1024*1024);
}

void double_array_init (double* p, size_t nelems)
{
    if (init_snan) array_init_snan(p, nelems);
}

void array_init_snan (double* p, size_t nelems)
{
    static_assert(sizeof(double) == sizeof(long long), "MemPool: sizeof double != sizeof long long");

    for (size_t i = 0; i < nelems; ++i) {
	long long *ll = (long long *) (p++);
	*ll = 0x7ff0000080000001LL;
    }
}

}
