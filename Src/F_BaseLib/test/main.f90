program main
  use BoxLib
  use parallel
  use multifab_module
  use bl_prof_module

  implicit none

  call boxlib_initialize()

  call bl_prof_initialize(on = .true.)

! call t_box_read
! call t_plotfile
! call t_ba
! call t_bx
! call t_cluster
! call t_mf_fabio
! call t_mf
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
! call t_ml_mf_read
! call t_bl_prof
! call t_ba_self_intersection
! call t_knapsack
! call t_gatherv
  call t_particle

 call layout_flush_copyassoc_cache ()

  if ( parallel_IOProcessor() ) then
     print*, 'MEMORY STATS AT END OF PROGRAM'
     print*, ' '
  end if
  call print(lmultifab_mem_stats(),   "   lmultifab")
  call print(multifab_mem_stats(),    "    multifab")
  call print(fab_mem_stats(),         "         fab")
  call print(boxarray_mem_stats(),    "    boxarray")
  call print(layout_mem_stats(),      "      layout")
  call print(boxassoc_mem_stats(),    "    boxassoc")
  call print(fgassoc_mem_stats(),     "     fgassoc")
  call print(syncassoc_mem_stats(),   "   syncassoc")
  call print(copyassoc_mem_stats(),   "   copyassoc")
  call print(fluxassoc_mem_stats(),   "   fluxassoc")

  call bl_prof_glean("bl_prof_res")
  call bl_prof_finalize()

  call boxlib_finalize()

end program main
