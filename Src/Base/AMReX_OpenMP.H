#ifndef AMREX_OPENMP_H_
#define AMREX_OPENMP_H_
#include <AMReX_Config.H>

#ifdef AMREX_USE_OMP
#include <omp.h>

namespace amrex::OpenMP {

    inline int get_num_threads () { return omp_get_num_threads(); }
    inline int get_max_threads () { return omp_get_max_threads(); }
    inline int get_thread_num  () { return omp_get_thread_num();  }
    inline int in_parallel     () { return omp_in_parallel();     }

}

#else

namespace amrex::OpenMP {

    constexpr int get_num_threads () { return 1; }
    constexpr int get_max_threads () { return 1; }
    constexpr int get_thread_num  () { return 0; }
    constexpr int in_parallel     () { return false; }

}

#endif

#endif
