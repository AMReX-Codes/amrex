program main

  implicit none

  logical ::  l
  integer, external :: omp_get_thread_num

  l = .true.
  if (l) then
     !$omp parallel
     print *, 'Thread #', omp_get_thread_num()
     !$omp end parallel
     ! ROSE will move the above line after "end if".
  end if

end program main
