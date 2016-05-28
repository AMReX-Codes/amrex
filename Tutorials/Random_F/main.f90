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
  type(bl_rng_binomial) :: bi_orig, bi_restore

  call boxlib_initialize()

  ! uniform real distribution: [0.0d0, 1.0d0), seed is 42.
  call bl_rng_build(ur_orig, 42, 0.0d0, 1.0d0)

  ! normal distribution with mean=0.d0, stddev=1.0d0, seed is 549
  call bl_rng_build(nm_orig, 549, 0.0d0, 1.0d0)

  ! poisson distribution with mean=100.d0, seed is 342
  call bl_rng_build(ps_orig, 342, 100.0d0)

  ! binomial distribution with t=160 (number of trials), p=0.6d0, seed is 8856
  call bl_rng_build(bi_orig, 8856, 160, 0.6d0)

  if (parallel_myproc() .eq. 0) then
     print *, 'uniform real, normal, poisson, binomial'
  end if
  do i = 1, 10
     r(1)  = bl_rng_get(ur_orig)
     r(2)  = bl_rng_get(nm_orig)
     ir(1) = bl_rng_get(ps_orig)
     ir(2) = bl_rng_get(bi_orig)

     if (parallel_myproc() .eq. 0) &
          print *, r(1), r(2), ir(1), ir(2)
  end do

  call bl_rng_save   (ur_orig   , "rng_state_uniform_real")
  call bl_rng_restore(ur_restore, "rng_state_uniform_real")

  call bl_rng_save   (nm_orig   , "rng_state_normal")
  call bl_rng_restore(nm_restore, "rng_state_normal")

  call bl_rng_save   (ps_orig   , "rng_state_poisson")
  call bl_rng_restore(ps_restore, "rng_state_poisson")

  call bl_rng_save   (bi_orig   , "rng_state_binomial")
  call bl_rng_restore(bi_restore, "rng_state_binomial")

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

  if (parallel_myproc() .eq. 0) then
     print *, 'binomial: original, restart, difference'
  end if
  do i = 1, 10
     ia = bl_rng_get(bi_orig)
     ib = bl_rng_get(bi_restore)
     if (parallel_myproc() .eq. 0) then
        print *, ia, ib, ia-ib
     else if (ia .ne. ib) then
        print *, "binomial error!!! ", ia, ib
     end if
  end do

  if (parallel_myproc() .eq. 0) then
     print *, 'rank, uniform real, normal, poisson, binomial'
  end if
  do i = 0, parallel_nprocs()
     if (i .eq. parallel_myproc()) then
        print *, parallel_myproc(), bl_rng_get(ur_orig), bl_rng_get(nm_orig), &
             bl_rng_get(ps_orig), bl_rng_get(bi_orig) 
        call flush(6)
     end if
     call parallel_barrier()
  end do

  if (parallel_myproc() .eq. 0) then
     print *, "we now change the change Poisson mean from 100.d0 to 30000.54d0"
     call bl_rng_change_distribution(ps_orig, 30000.54d0)
     print *, "Poisson old and new"
     do i = 1, 5
        print *, bl_rng_get(ps_restore), bl_rng_get(ps_orig)
     end do
  end if
  
  if (parallel_myproc() .eq. 0) then
     print *, "we now change the change Binomial mean from 160, 0.6d0 to 1000, 0.3d0"
     call bl_rng_change_distribution(bi_orig, 1000, 0.3d0)
     print *, "Binomial old and new"
     do i = 1, 5
        print *, bl_rng_get(bi_restore), bl_rng_get(bi_orig)
     end do
  end if
  
  call bl_rng_destroy(ur_orig)
  call bl_rng_destroy(ur_restore)
  call bl_rng_destroy(nm_orig)
  call bl_rng_destroy(nm_restore)
  call bl_rng_destroy(ps_orig)
  call bl_rng_destroy(ps_restore)
  call bl_rng_destroy(bi_orig)
  call bl_rng_destroy(bi_restore)

  call boxlib_finalize()

end program main
