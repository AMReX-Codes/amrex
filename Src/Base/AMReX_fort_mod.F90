module amrex_fort_module

  use iso_c_binding, only : c_float, c_double, c_size_t, c_ptr

  implicit none

  integer, parameter ::    bl_spacedim = AMREX_SPACEDIM
  integer, parameter :: amrex_spacedim = AMREX_SPACEDIM

#ifdef BL_USE_FLOAT
  integer, parameter :: amrex_real = c_float
  ! We could/should use Fortran 2008 c_sizeof here.
  integer (kind=c_size_t), parameter :: amrex_real_size = 4_c_size_t
#else
  integer, parameter :: amrex_real = c_double
  ! We could/should use Fortran 2008 c_sizeof here.
  integer (kind=c_size_t), parameter :: amrex_real_size = 8_c_size_t
#endif

#ifdef BL_SINGLE_PRECISION_PARTICLES
  integer, parameter :: amrex_particle_real = c_float
#else
  integer, parameter :: amrex_particle_real = c_double
#endif

  interface
     function amrex_malloc (s) bind(c,name='amrex_malloc')
       import
       integer(c_size_t), intent(in), value :: s
       type(c_ptr) :: amrex_malloc
     end function amrex_malloc

     subroutine amrex_free (p) bind(c,name='amrex_free')
       import
       type(c_ptr), value :: p
     end subroutine amrex_free
  end interface

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



  AMREX_DEVICE subroutine get_loop_bounds(blo, bhi, lo, hi) bind(c, name='get_loop_bounds')

    implicit none

    integer, intent(in   ) :: lo(3), hi(3)
    integer, intent(inout) :: blo(3), bhi(3)

#ifdef AMREX_USE_CUDA
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



  AMREX_DEVICE subroutine amrex_add(x, y)

    implicit none

    ! Add y to x. Do it safely if we're on the GPU.

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    real(amrex_real) :: t

#ifdef AMREX_USE_CUDA
    t = atomicAdd(x, y)
#else
    x = x + y
#endif

  end subroutine amrex_add



  AMREX_DEVICE subroutine amrex_subtract(x, y)

    implicit none

    ! Subtract y from x. Do it safely if we're on the GPU.

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    real(amrex_real) :: t

#ifdef AMREX_USE_CUDA
    t = atomicSub(x, y)
#else
    x = x - y
#endif

  end subroutine amrex_subtract



  AMREX_DEVICE subroutine amrex_max(x, y)

    implicit none

    ! Set in x the maximum of x and y. Do it safely if we're on the GPU.

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    real(amrex_real) :: t

#ifdef AMREX_USE_CUDA
    t = atomicMax(x, y)
#else
    x = max(x, y)
#endif

  end subroutine amrex_max



  AMREX_DEVICE subroutine amrex_min(x, y)

    implicit none

    ! Set in x the minimum of x and y. Do it safely if we're on the GPU.

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    real(amrex_real) :: t

#ifdef AMREX_USE_CUDA
    t = atomicMin(x, y)
#else
    x = min(x, y)
#endif

  end subroutine amrex_min

end module amrex_fort_module
