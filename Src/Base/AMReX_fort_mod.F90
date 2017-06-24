module amrex_fort_module

  use iso_c_binding, only : c_float, c_double

  implicit none

  integer, parameter ::    bl_spacedim = BL_SPACEDIM
  integer, parameter :: amrex_spacedim = BL_SPACEDIM

#ifdef BL_USE_FLOAT
  integer, parameter :: amrex_real = c_float
#else
  integer, parameter :: amrex_real = c_double
#endif

#ifdef BL_SINGLE_PRECISION_PARTICLES
  integer, parameter :: amrex_particle_real = c_float
#else
  integer, parameter :: amrex_particle_real = c_double
#endif

contains

  function amrex_coarsen_intvect (n, iv, rr) result(civ)
    integer, intent(in) :: n, rr
    integer, intent(in) :: iv(n)
    integer :: civ(n)
    integer :: i
    do i = 1, n
       if (iv(i) .lt. 0) then
          civ(i) = -abs(iv(i)+1)/rr - 1
       else
          civ(i) = iv(i)/rr
       end if
    end do
  end function amrex_coarsen_intvect

#ifdef CUDA
  attributes(device) &
#endif
  subroutine get_loop_bounds(blo, bhi, lo, hi)

    implicit none

    integer, intent(in   ) :: lo(3), hi(3)
    integer, intent(inout) :: blo(3), bhi(3)

#ifdef CUDA
    ! Get our spatial index based on the CUDA thread index

    blo(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    blo(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    blo(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    ! If we have more threads than zones, set hi < lo so that the
    ! loop iteration gets skipped.

    if (blo(1) .gt. hi(1) .or. blo(2) .gt. hi(2) .or. blo(3) .gt. hi(3)) then
       bhi = blo - 1
    else
       bhi = blo
    endif
#else
    blo = lo
    bhi = hi
#endif

  end subroutine get_loop_bounds

end module amrex_fort_module
