#ifndef AMREX_OPENMP_H_
#define AMREX_OPENMP_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_OpenMP.H
 * \brief OpenMP utility wrappers and thread-count query functions.
 */

#ifdef AMREX_USE_OMP
#include <AMReX_Extension.H>
#include <omp.h>

namespace amrex::OpenMP {

    //! Return the number of threads in the current parallel region.
    inline int get_num_threads () { return omp_get_num_threads(); }
    //! Return the maximum number of threads that may be used in a parallel region.
    inline int get_max_threads () { return omp_get_max_threads(); }
    //! Return the thread number of the calling thread within the current parallel region.
    inline int get_thread_num  () { return omp_get_thread_num();  }
    //! Return non-zero if called from within an active parallel region.
    inline int in_parallel     () { return omp_in_parallel();     }
    //! Set the number of threads to use for subsequent parallel regions.
    inline void set_num_threads (int num) { omp_set_num_threads(num); }

    //! Initialize the OpenMP subsystem. Used internally by AMReX.
    void Initialize ();
    //! Finalize the OpenMP subsystem. Used internally by AMReX.
    void Finalize ();

#if defined(_WIN32)
    /// \cond DOXYGEN_IGNORE
    void** get_lock_impl (int ilock);
    /// \endcond

    //! Return the OpenMP lock at index \p ilock.
    inline omp_lock_t* get_lock (int ilock) {
        return (omp_lock_t*)(*(get_lock_impl(ilock)));
    }
#else
    //! Return the OpenMP lock at index \p ilock.
    omp_lock_t* get_lock (int ilock);
#endif
}

#else // AMREX_USE_OMP

namespace amrex::OpenMP {

    //! Return the number of threads in the current parallel region (always 1 without OpenMP).
    constexpr int get_num_threads () { return 1; }
    //! Return the maximum number of threads (always 1 without OpenMP).
    constexpr int get_max_threads () { return 1; }
    //! Return the thread number of the calling thread (always 0 without OpenMP).
    constexpr int get_thread_num  () { return 0; }
    //! Return whether in a parallel region (always false without OpenMP).
    constexpr int in_parallel     () { return false; }
    //! No-op: set the number of threads (no-op without OpenMP).
    constexpr void set_num_threads (int) { /* nothing */ }
}

#endif // AMREX_USE_OMP

namespace amrex {
    /// \ingroup amrex_utilities
    //! Return the number of unique physical CPU cores available to this process.
    int numUniquePhysicalCores();
}

#endif
