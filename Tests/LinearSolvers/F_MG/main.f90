program main
  use BoxLib
  use bl_prof_module
  use multifab_module
  use layout_module

  implicit none

  call boxlib_initialize()
  call bl_prof_initialize(on = .true.)

!  call layout_set_verbosity(1)
  call layout_set_copyassoc_max(25)

  !  call multifab_set_behavior(fb_async = .true., fb_fancy = .true.)

  ! call t_cc_multigrid()
  ! call t_nodal_multigrid()

  call wrapper()
  call layout_flush_copyassoc_cache ()

  ! call t_stencil
  ! call t_stencil_map
  ! call t_smoother

  call bl_prof_glean("bl_prof_res")
  call bl_prof_finalize()
  call boxlib_finalize()

end program main
