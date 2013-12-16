program main
  use BoxLib
  use bl_prof_module
  use layout_module
  use multifab_module

  implicit none

  call boxlib_initialize()
  call bl_prof_initialize(on = .true.)

  call layout_set_copyassoc_max(25)

  call wrapper()
  call layout_flush_copyassoc_cache ()

!  call print(multifab_mem_stats(),  " multifab at end")
!  call print(imultifab_mem_stats(), "imultifab at end")
!  call print(fab_mem_stats(),       "      fab at end")
!  call print(ifab_mem_stats(),      "     ifab at end")
!  call print(boxarray_mem_stats(),  " boxarray at end")
!  call print(boxassoc_mem_stats(),  " boxassoc at end")
!  call print(layout_mem_stats(),    "   layout at end")

  call bl_prof_glean("bl_prof_res")
  call bl_prof_finalize()
  call boxlib_finalize()

end program main
