program main

  use BoxLib
  use parallel
  use bl_random_module
  use iso_c_binding

  implicit none
  
  integer :: i
  double precision :: ra, rb, r(2)
  integer :: ia, ib, ir(2)
  type(bl_rng_uniform_real) :: ur_orig, ur_restore
  type(bl_rng_normal) :: nm_orig, nm_restore
  type(bl_rng_poisson) :: ps_orig, ps_restore

  call boxlib_initialize()

  ! uniform real distribution: [0.0d0, 1.0d0), seed is 42.
  call bl_rng_build(ur_orig, 42, 0.0d0, 1.0d0)

  ! normal distribution with mean=0.d0, stddev=1.0d0
  call bl_rng_build(nm_orig, 549, 0.0d0, 1.0d0)

  ! poisson distribution with mean=10.d0
  call bl_rng_build(ps_orig, 342, 10.0d0)

  if (parallel_myproc() .eq. 0) then
     print *, 'uniform real, normal, poisson'
  end if
  do i = 1, 10
     r(1)  = bl_rng_get(ur_orig)
     r(2)  = bl_rng_get(nm_orig)
     ir(1) = bl_rng_get(ps_orig)

     if (parallel_myproc() .eq. 0) &
          print *, r(1), r(2), ir(1)
  end do

  call bl_rng_save   (ur_orig   , "rng_state_uniform_real")
  call bl_rng_restore(ur_restore, "rng_state_uniform_real")

  call bl_rng_save   (nm_orig   , "rng_state_normal")
  call bl_rng_restore(nm_restore, "rng_state_normal")

  call bl_rng_save   (ps_orig   , "rng_state_poisson")
  call bl_rng_restore(ps_restore, "rng_state_poisson")


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

  if (parallel_myproc() .eq. 0) then
     print *, 'normal: original, restart, difference'
  end if
  do i = 1, 10
     ra = bl_rng_get(nm_orig)
     rb = bl_rng_get(nm_restore)
     if (parallel_myproc() .eq. 0) then
        print *, ra, rb, ra-rb
     else if (ra .ne. rb) then
        print *, "normal error!!! ", ra, rb
     end if
  end do

  if (parallel_myproc() .eq. 0) then
     print *, 'poisson: original, restart, difference'
  end if
  do i = 1, 10
     ia = bl_rng_get(ps_orig)
     ib = bl_rng_get(ps_restore)
     if (parallel_myproc() .eq. 0) then
        print *, ia, ib, ia-ib
     else if (ia .ne. ib) then
        print *, "poisson error!!! ", ia, ib
     end if
  end do
  
  call bl_rng_destroy(ur_orig)
  call bl_rng_destroy(ur_restore)
  call bl_rng_destroy(nm_orig)
  call bl_rng_destroy(nm_restore)
  call bl_rng_destroy(ps_orig)
  call bl_rng_destroy(ps_restore)

  call boxlib_finalize()

end program main
