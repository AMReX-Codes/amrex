#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>

#include <CArena.H>
#include <PArray.H>
#include <MemPool.H>

namespace
{
  static PArray<CArena> memory_pool;
}

void init_memmory_pool()
{
#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
#else
  int nthreads = 1;
#endif
  memory_pool.resize(nthreads, PArrayManage);
  for (int i=0; i<nthreads; ++i) {
    memory_pool.set(i, new CArena());
  }
}

extern "C"
void* bl_allocate_c (size_t nbytes)
{
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  return memory_pool[tid].alloc(nbytes);
}

extern "C"
void bl_deallocate_c (void* p) 
{
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  memory_pool[tid].free(p);
}
