module omp_module

  implicit none

  integer, external :: omp_get_num_threads
  integer, external :: omp_get_max_threads
  integer, external :: omp_get_thread_num
  integer, external :: omp_get_num_procs
  external omp_set_num_threads
  external omp_set_dynamic
  logical, external :: omp_get_dynamic
  logical, external :: omp_in_parallel
  external omp_set_nested
  logical, external :: omp_get_nested
  external omp_init_lock
  external omp_destroy_lock
  external omp_set_lock
  external omp_unset_lock
  logical, external :: omp_test_lock

end module omp_module
