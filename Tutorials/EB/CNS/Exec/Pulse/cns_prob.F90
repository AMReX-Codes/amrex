module probdata_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  real(rt), save :: rpulse = 0.5d0
  real(rt), save :: rho0   = 1.2d-3
  real(rt), save :: drho0  = 1.2d-4
  real(rt), save :: p0     = 1.01325d6
end module probdata_module


subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
  use amrex_fort_module, only : rt => amrex_real
  use amrex_parmparse_module
  use probdata_module
  implicit none
  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(*), probhi(*)
  type(amrex_parmparse) :: pp
  call amrex_parmparse_build(pp,"prob")
  call pp%query("rpulse",rpulse)
  call pp%query("rho",rho0)
  call pp%query("drho",drho0)
  call pp%query("p0",p0)
  call amrex_parmparse_destroy(pp)
end subroutine amrex_probinit


subroutine cns_initdata(level, time, lo, hi, u, ulo, uhi, dx, prob_lo) bind(C, name="cns_initdata")
  use amrex_fort_module, only : rt => amrex_real
  use cns_physics_module, only : gamma, cv
  use cns_module, only : center, nvar, urho, umx, umy, umz, ueden, ueint, utemp
  use probdata_module, only : rpulse, rho0, drho0, p0
  implicit none
  integer, intent(in) :: level, lo(3), hi(3), ulo(3), uhi(3)
  real(rt), intent(in) :: time
  real(rt), intent(inout) :: u(ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3),nvar)
  real(rt), intent(in) :: dx(3), prob_lo(3)
  
  integer :: i,j,k
  real(rt) :: x,y,z,r, Pt
  real(rt), parameter :: Pi = 4.d0*atan(1.d0)

  do k = lo(3), hi(3)
     z = prob_lo(3) + (k+0.5d0)*dx(3) - center(3)
     do j = lo(2), hi(2)
        y = prob_lo(2) + (j+0.5d0)*dx(2) - center(2)
        do i = lo(1), hi(1)
           x = prob_lo(1) + (i+0.5d0)*dx(1) - center(1)
           r = sqrt(x*x + y*y + z*z)

           if (r .gt. rpulse) then
              u(i,j,k,urho) = rho0
           else
              u(i,j,k,urho) = rho0 + drho0*exp(-16.d0*r*r)*(cos(Pi*r))**6
           end if

           Pt = p0*(u(i,j,k,urho)/rho0)**gamma
           u(i,j,k,ueint) = Pt / (gamma-1.d0)

           u(i,j,k,umx:umz) = 0.d0
           u(i,j,k,ueden) = u(i,j,k,ueint)
           u(i,j,k,utemp) = u(i,j,k,ueint)/(u(i,j,k,urho)*cv)

        end do
     end do
  end do

end subroutine cns_initdata
