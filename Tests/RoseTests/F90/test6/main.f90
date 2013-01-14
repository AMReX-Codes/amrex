program main

  use omp_module

  implicit none

  print *, "Max threads: ", omp_get_max_threads(), " threads"

end program main
