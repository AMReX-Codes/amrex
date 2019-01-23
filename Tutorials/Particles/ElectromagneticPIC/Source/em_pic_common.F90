module constants
  use amrex_fort_module, only : amrex_real

  real(amrex_real), parameter :: clight  = 2.99792458d8
  real(amrex_real), parameter :: epsilon_0 = 8.85418782d-12
  real(amrex_real), parameter :: electron_charge = 1.60217662d-19
  real(amrex_real), parameter :: electron_mass = 9.10938356d-31

end module

module em_particle_module
  use amrex_fort_module, only: amrex_real, amrex_particle_real
  use iso_c_binding ,    only: c_int

  implicit none
  private

  public  particle_t

  type, bind(C)  :: particle_t
     real(amrex_particle_real) :: pos(3)     !< Position
     integer(c_int)            :: id         !< Particle id
     integer(c_int)            :: cpu        !< Particle cpu
  end type particle_t

end module em_particle_module

subroutine check_langmuir_solution(boxlo, boxhi, testlo, testhi, jx, jxlo, jxhi, &
     time, max_error) bind(c,name='check_langmuir_solution')

  use amrex_fort_module, only : amrex_real
  use constants

  implicit none

  integer,          intent(in)    :: boxlo(3),  boxhi(3)
  integer,          intent(in)    :: testlo(3), testhi(3)
  integer,          intent(in)    :: jxlo(3),   jxhi(3)
  real(amrex_real), intent(in)    :: jx(jxlo(1):jxhi(1),jxlo(2):jxhi(2),jxlo(3):jxhi(3))
  real(amrex_real), intent(in), value :: time
  real(amrex_real), intent(inout) :: max_error

  integer :: j,k,l
  integer :: lo(3), hi(3)

  real(amrex_real) error
  real(amrex_real) exact

  real(amrex_real), parameter :: u  = 0.01
  real(amrex_real), parameter :: n0 = 1.d25
  real(amrex_real) wp

  wp = sqrt(n0*electron_charge**2/(electron_mass*epsilon_0))

  error = 0.0
  exact = -n0*electron_charge*clight*u*cos(wp*time)

  lo = max(boxlo, testlo)
  hi = min(boxhi, testhi)

  do l       = lo(3), hi(3)
     do k    = lo(2), hi(2)
        do j = lo(1), hi(1)
           error = max(error, abs(jx(j,k,l) - exact) / abs(exact))
        end do
     end do
  end do

  max_error = error

end subroutine check_langmuir_solution

