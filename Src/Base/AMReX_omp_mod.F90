#ifdef BL_USE_OMP

module amrex_omp_module

  implicit none

#ifdef _OPENMP
  integer, parameter :: amrex_omp_support = (_OPENMP)
#else
  integer, parameter :: amrex_omp_support = 0 ! Should this be allowed??
#endif

  integer, external :: omp_get_num_threads
  integer, external :: omp_get_max_threads
  integer, external :: omp_get_thread_num
  logical, external :: omp_in_parallel

end module amrex_omp_module

#else

module amrex_omp_module

  implicit none

  integer, parameter :: amrex_omp_support = 0 ! indicates no support

contains

  integer function omp_get_num_threads()
    omp_get_num_threads = 1
  end function omp_get_num_threads

  integer function omp_get_max_threads()
    omp_get_max_threads = 1
  end function omp_get_max_threads

  integer function omp_get_thread_num()
    omp_get_thread_num = 0
  end function omp_get_thread_num

  logical function omp_in_parallel()
    omp_in_parallel = .false.
  end function omp_in_parallel

end module amrex_omp_module

#endif

