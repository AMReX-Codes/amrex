module uniform01_module

  use amrex_fort_module, only : rt => amrex_real
  use bl_random_module

  implicit none
  private

  public :: get_uniform01

  type(bl_rng_engine), save :: eng_uniform01
  type(bl_rng_uniform_real), save :: dist_uniform01
  !$OMP THREADPRIVATE (eng_uniform01, dist_uniform01)
  
  contains

! :::
! ::: ----------------------------------------------------------------
! :::

  subroutine init_uniform01_rng() &
       bind(c,name='init_uniform01_rng')

    use iso_c_binding

    ! Right now, bl_random_module is NOT thread-safe:
    ! every thread on the same processor gets the same random numbers.

    ! seed 0 picks a random root seed
    !$OMP parallel
    call bl_rng_build_engine(eng_uniform01, 0)

    ! uniform real distribution: [0.0d0, 1.0d0)
    call bl_rng_build_distro(dist_uniform01, 0.0d0, 1.0d0)
    !$OMP end parallel
    
  end subroutine init_uniform01_rng

! :::
! ::: ----------------------------------------------------------------
! :::

  function get_uniform01() result(r)
    real(rt) :: r
    r = bl_rng_get(dist_uniform01, eng_uniform01)
  end function get_uniform01

end module uniform01_module
