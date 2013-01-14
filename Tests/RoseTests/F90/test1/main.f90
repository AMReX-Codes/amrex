program main

  implicit none

  logical :: pmask(3)

  pmask = .false.

  ! ROSE translates .eqv. in the following line to .eq.
  if ( all(pmask .eqv. .false.) ) then
     print *, 'Hello world!'
  end if

  ! ROSE translates .eqv. in the following line to .eq.
  if ( any(pmask .eqv. .true.) ) then
     print *, 'How did this happen?'
  end if

end program main
