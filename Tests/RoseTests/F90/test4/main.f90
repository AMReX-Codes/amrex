program main

  implicit none

  ! ROSE translates this to 
  ! REAL(kind=8) :: x = __builtin_huge_valf()
  real(kind=8) :: x = 1.e50_8

end program main
