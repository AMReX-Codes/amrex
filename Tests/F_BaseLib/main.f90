program main
  use BoxLib
  implicit none

  call boxlib_initialize()

! call t_pingpong
! call t_plotfile
! call t_ba
call t_cluster
! call t_mf_fabio
! call t_mf
! call t_fabio
! call t_domain
! call t_box_chop
! call t_box_mod
! call t_boxassoc
! call t_mt_random_numbers

  call boxlib_finalize()

end program main
