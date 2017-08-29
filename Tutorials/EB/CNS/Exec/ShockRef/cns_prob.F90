module probdata_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  real(rt), save :: p0   = 1.95d4
  real(rt), save :: p1   = 7.42d6
  real(rt), save :: rho0 = 3.29d-5
  real(rt), save :: rho1 = 3.61d-4
  real(rt), save :: v0   = 0.d0
  real(rt), save :: v1   = 0.d0
  real(rt), save :: x1   = 4.2
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
  call pp%query("p0", p0)
  call pp%query("p1", p1)
  call pp%query("rho0", rho0)
  call pp%query("rho1", rho1)
  call pp%query("v0", v0)
  call pp%query("v1", v1)
  call pp%query("x1", x1)
  call amrex_parmparse_destroy(pp)
end subroutine amrex_probinit


subroutine cns_initdata(level, time, lo, hi, u, ulo, uhi, dx, prob_lo) bind(C, name="cns_initdata")
  use amrex_fort_module, only : rt => amrex_real
  use cns_physics_module, only : gamma, cv
  use cns_module, only : center, nvar, urho, umx, umy, umz, ueden, ueint, utemp
  use probdata_module, only : p0,p1,rho0,rho1,v0,v1,x1
  implicit none
  integer, intent(in) :: level, lo(3), hi(3), ulo(3), uhi(3)
  real(rt), intent(in) :: time
  real(rt), intent(inout) :: u(ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3),nvar)
  real(rt), intent(in) :: dx(3), prob_lo(3)
  
  integer :: i,j,k
  real(rt) :: x

  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)

           x = (i+0.5d0)*dx(1)

           if (x .gt. x1) then
              u(i,j,k,urho) = rho0
              u(i,j,k,ueint) = p0 / (gamma-1.d0)
              u(i,j,k,umx) = u(i,j,k,urho)*v0
           else
              u(i,j,k,urho) = rho1
              u(i,j,k,ueint) = p1 / (gamma-1.d0)
              u(i,j,k,umx) = u(i,j,k,urho)*v1
           end if

           u(i,j,k,umy:umz) = 0.d0
           u(i,j,k,ueden) = u(i,j,k,ueint) + 0.5d0*u(i,j,k,umx)**2/u(i,j,k,urho)
           u(i,j,k,utemp) = u(i,j,k,ueint)/(u(i,j,k,urho)*cv)

        end do
     end do
  end do

end subroutine cns_initdata
