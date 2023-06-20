#include <AMReX_Config.H>

module amrex_fort_module

  use iso_c_binding, only : c_char, c_short, c_int, c_long, c_long_long, c_float, c_double, c_size_t, c_ptr

  implicit none

  integer, parameter ::    bl_spacedim = AMREX_SPACEDIM
  integer, parameter :: amrex_spacedim = AMREX_SPACEDIM

#ifdef AMREX_USE_FLOAT
  integer, parameter :: amrex_real = c_float
  ! We could/should use Fortran 2008 c_sizeof here.
  integer (kind=c_size_t), parameter :: amrex_real_size = 4_c_size_t
#else
  integer, parameter :: amrex_real = c_double
  ! We could/should use Fortran 2008 c_sizeof here.
  integer (kind=c_size_t), parameter :: amrex_real_size = 8_c_size_t
#endif

#ifdef AMREX_SINGLE_PRECISION_PARTICLES
  integer, parameter :: amrex_particle_real = c_float
#else
  integer, parameter :: amrex_particle_real = c_double
#endif

#ifdef _WIN32
  integer, parameter :: amrex_long = c_long_long
#else
  integer, parameter :: amrex_long = c_long;
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

     function amrex_random () bind(c,name='amrex_random')
       import
       real(c_double) :: amrex_random
     end function amrex_random

     function amrex_random_int (n) bind(c,name='amrex_random_int')
       import
       integer(amrex_long), intent(in), value :: n
       integer(amrex_long) :: amrex_random_int
     end function amrex_random_int
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

end module amrex_fort_module
