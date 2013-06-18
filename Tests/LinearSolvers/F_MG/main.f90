program main
  use BoxLib
  use bl_prof_module
  use layout_module

  implicit none

  call boxlib_initialize()
  call bl_prof_initialize(on = .true.)

  call layout_set_copyassoc_max(25)

  call wrapper()
  call layout_flush_copyassoc_cache ()

  call bl_prof_glean("bl_prof_res")
  call bl_prof_finalize()
  call boxlib_finalize()

end program main
