program main

  use BoxLib
  use parallel
  use bl_random_module
  use iso_c_binding

  implicit none
  
  integer :: i
  double precision :: ra, rb, r(2)
  integer :: ia, ib, ir(2)
  type(bl_rng_engine) :: eng1, eng2, eng1_r, eng2_r
  type(bl_rng_uniform_real) :: ur, ur_r
  type(bl_rng_normal) :: nm, nm_r
  type(bl_rng_poisson) :: ps, ps_r
  type(bl_rng_binomial) :: bi, bi_r

  call boxlib_initialize()

  call bl_rng_build_engine(eng1, 42)
  call bl_rng_build_engine(eng2, 0)

  ! uniform real distribution: [0.0d0, 1.0d0)
  call bl_rng_build_distro(ur, 0.0d0, 1.0d0)

  ! normal distribution with mean=0.d0, stddev=1.0d0
  call bl_rng_build_distro(nm, 0.0d0, 1.0d0)

  ! poisson distribution with mean=100.d0
  call bl_rng_build_distro(ps, 100.0d0)

  ! binomial distribution with t=160 (number of trials), p=0.6d0
  call bl_rng_build_distro(bi, 160, 0.6d0)

  if (parallel_myproc() .eq. 0) then
     print *, 'uniform real, normal, poisson, binomial'
  end if
  do i = 1, 10
     r(1)  = bl_rng_get(ur, eng1)
     r(2)  = bl_rng_get(nm, eng2)
     ir(1) = bl_rng_get(ps, eng2)
     ir(2) = bl_rng_get(bi, eng2)

     if (parallel_myproc() .eq. 0) &
          print *, r(1), r(2), ir(1), ir(2)
  end do

  call bl_rng_save_engine   (eng1  , "eng1")
  call bl_rng_restore_engine(eng1_r, "eng1")

  call bl_rng_save_engine   (eng2  , "eng2")
  call bl_rng_restore_engine(eng2_r, "eng2")

  call bl_rng_save_distro   (nm  , "normal")
  call bl_rng_restore_distro(nm_r, "normal")

  call bl_rng_save_distro   (ur  , "uniform")
  call bl_rng_restore_distro(ur_r, "uniform")

  call bl_rng_save_distro   (ps  , "poisson")
  call bl_rng_restore_distro(ps_r, "poisson")

  call bl_rng_save_distro   (bi  , "binomial")
  call bl_rng_restore_distro(bi_r, "binomial")

  if (parallel_myproc() .eq. 0) then
     print *, 'uniform real: original, restart, difference'
  end if
  do i = 1, 10
     ra = bl_rng_get(ur  , eng2)
     rb = bl_rng_get(ur_r, eng2_r)
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
     ra = bl_rng_get(nm  , eng2)
     rb = bl_rng_get(nm_r, eng2_r)
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
     ia = bl_rng_get(ps  , eng1)
     ib = bl_rng_get(ps_r, eng1_r)
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
     ia = bl_rng_get(bi  , eng2)
     ib = bl_rng_get(bi_r, eng2_r)
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
        print *, parallel_myproc(), bl_rng_get(ur,eng1), bl_rng_get(nm,eng2), &
             bl_rng_get(ps,eng2), bl_rng_get(bi,eng1) 
        flush(6)
     end if
     call parallel_barrier()
  end do

  if (parallel_myproc() .eq. 0) then
     print *, "we now change the change Poisson mean from 100.d0 to 30000.54d0"
     call bl_rng_destroy_distro(ps)
     call bl_rng_build_distro(ps, 30000.54d0)
     print *, "Poisson old and new"
     do i = 1, 5
        print *, bl_rng_get(ps_r, eng1), bl_rng_get(ps, eng2)
     end do
  end if
  
  if (parallel_myproc() .eq. 0) then
     print *, "we now change the change Binomial mean from 160, 0.6d0 to 1000, 0.3d0"
     call bl_rng_destroy_distro(bi)
     call bl_rng_build_distro(bi, 1000, 0.3d0)
     print *, "Binomial old and new"
     do i = 1, 5
        print *, bl_rng_get(bi_r, eng1), bl_rng_get(bi, eng2)
     end do
  end if
  
  call bl_rng_destroy_distro(ur)
  call bl_rng_destroy_distro(nm)
  call bl_rng_destroy_distro(ps)
  call bl_rng_destroy_distro(bi)
  call bl_rng_destroy_distro(ur_r)
  call bl_rng_destroy_distro(nm_r)
  call bl_rng_destroy_distro(ps_r)
  call bl_rng_destroy_distro(bi_r)

  call bl_rng_destroy_engine(eng1)
  call bl_rng_destroy_engine(eng2)
  call bl_rng_destroy_engine(eng1_r)
  call bl_rng_destroy_engine(eng2_r)

  call boxlib_finalize()

end program main
