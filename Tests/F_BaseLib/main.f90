program main
  use BoxLib
  use bl_prof_module
  implicit none
  external log2
  integer log2, i

  call boxlib_initialize()
  call bl_prof_initialize(on = .true.)

! call t_box_read()
! call t_pingpong
! call t_plotfile
! call t_ba
! call t_bx
! call t_cluster
! call t_mf_fabio
! call t_mf
! call t_fabio
! call t_domain
! call t_box_chop
! call t_box_mod
! call t_boxassoc
! call t_mt_random_numbers
! call t_box_conn
! call t_boxarray
! call t_nodal_mf_fabio
! call t_timer
! call t_knap
! call t_ml_mf_read()
! call t_bl_prof()
! call t_ba_self_intersection
  call t_knapsack
  call bl_prof_glean("bl_prof_res")
  call bl_prof_finalize()
  call boxlib_finalize()

end program main
