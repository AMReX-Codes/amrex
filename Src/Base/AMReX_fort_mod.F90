module amrex_fort_module

  use iso_c_binding, only : c_float, c_double, c_size_t

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
