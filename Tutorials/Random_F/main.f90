program main

  use BoxLib
  use parallel
  use bl_random_module
  use iso_c_binding

  implicit none
  
  integer :: i
  double precision :: ra, rb
  type(bl_rng_uniform_real) :: ur_orig, ur_restore
  type(bl_rng_normal) :: nm_orig, nm_restore

  call boxlib_initialize()

  ! uniform real distribution: [0.0d0, 1.0d0), seed is 42.
  call bl_rng_build(ur_orig, 42, 0.0d0, 1.0d0)

  if (parallel_myproc() .eq. 0) then
     do i = 1, 10
        print *, "uniform real: ", bl_rng_get(ur_orig)
     end do
  end if

  call bl_rng_save(ur_orig, "rng_state_uniform_real")

  call bl_rng_restore(ur_restore, "rng_state_uniform_real")

  if (parallel_myproc() .eq. 0) then
     print *, 'uniform real: original, restart, difference'
  end if
  do i = 1, 10
     ra = bl_rng_get(ur_orig)
     rb = bl_rng_get(ur_restore)
     if (parallel_myproc() .eq. 0) then
        print *, ra, rb, ra-rb
     else if (ra .ne. rb) then
        print *, "uniform real error!!! ", ra, rb
     end if
  end do

  ! normal distribution with mean=0.d0, stddev=1.0d0
  call bl_rng_build(nm_orig, 549, 0.0d0, 1.0d0)

  if (parallel_myproc() .eq. 0) then
     do i = 1, 10
        print *, "normal: ", bl_rng_get(nm_orig)
     end do
  end if

  call bl_rng_save(nm_orig, "rng_state_normal")

  call bl_rng_restore(nm_restore, "rng_state_normal")

  if (parallel_myproc() .eq. 0) then
     print *, 'normal: original, restart, difference'
  end if
  do i = 1, 10
     ra = bl_rng_get(nm_orig)
     rb = bl_rng_get(nm_restore)
     if (parallel_myproc() .eq. 0) then
        print *, ra, rb, ra-rb
     else if (ra .ne. rb) then
        print *, "nomal error!!! ", ra, rb
     end if
  end do
  
  call bl_rng_destroy(ur_orig)
  call bl_rng_destroy(ur_restore)
  call bl_rng_destroy(nm_orig)
  call bl_rng_destroy(nm_restore)

  call boxlib_finalize()

end program main
