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



  subroutine get_loop_bounds(blo, bhi, lo, hi)

    implicit none

    integer, intent(in   ) :: lo(3), hi(3)
    integer, intent(inout) :: blo(3), bhi(3)

    blo = lo
    bhi = hi

  end subroutine get_loop_bounds



  subroutine amrex_add(x, y)

    implicit none

    ! Add y to x.

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    x = x + y

  end subroutine amrex_add



  subroutine amrex_subtract(x, y)

    implicit none

    ! Subtract y from x.

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    x = x - y

  end subroutine amrex_subtract



  subroutine amrex_max(x, y)

    implicit none

    ! Set in x the maximum of x and y.

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    x = max(x, y)

  end subroutine amrex_max



  subroutine amrex_min(x, y)

    implicit none

    ! Set in x the minimum of x and y.

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    x = min(x, y)

  end subroutine amrex_min

end module amrex_fort_module
