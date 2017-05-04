module amrex_particle_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real

  implicit none

  private

  public :: amrex_particle_set_position, amrex_particle_get_position

contains

  subroutine amrex_particle_set_position (particles, ns, np, x) &
       bind(c,name='amrex_particle_set_position')
    integer(c_int)  , intent(in   ), value :: ns, np
    real(amrex_real), intent(inout)        :: particles(ns,np)
    real(amrex_real), intent(in   )        :: x(np)

    integer :: i

    do i = 1, np
       particles(1,i) = x(i)
    end do
  end subroutine amrex_particle_set_position

  subroutine amrex_particle_get_position (particles, ns, np, x) &
       bind(c,name='amrex_particle_get_position')
    integer(c_int)  , intent(in   ), value :: ns, np
    real(amrex_real), intent(in   )        :: particles(ns,np)
    real(amrex_real), intent(  out)        :: x(np)

    integer :: i

    do i = 1, np
       x(i) = particles(1,i)
    end do
  end subroutine amrex_particle_get_position

end module amrex_particle_module
