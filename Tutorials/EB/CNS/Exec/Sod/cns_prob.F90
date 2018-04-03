module probdata_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  real(rt), save :: p_l   = 1.0d0
  real(rt), save :: p_r   = 0.1d0
  real(rt), save :: rho_l = 1.0d0
  real(rt), save :: rho_r = 0.125d0
  real(rt), save :: u_l   = 0.d0
  real(rt), save :: u_r   = 0.d0
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
  call pp%query("p_l",p_l)
  call pp%query("p_r",p_r)
  call pp%query("rho_l",rho_l)
  call pp%query("rho_r",rho_r)
  call pp%query("u_l",u_l)
  call pp%query("u_r",u_r)
  call amrex_parmparse_destroy(pp)
end subroutine amrex_probinit


subroutine cns_initdata(level, time, lo, hi, u, ulo, uhi, dx, prob_lo) bind(C, name="cns_initdata")
  use amrex_fort_module, only : rt => amrex_real
  use cns_physics_module, only : gamma, cv
  use cns_module, only : nvar, urho, umx, umy, umz, ueden, ueint, utemp
  use probdata_module
  implicit none
  integer, intent(in) :: level, lo(3), hi(3), ulo(3), uhi(3)
  real(rt), intent(in) :: time
  real(rt), intent(inout) :: u(ulo(1):uhi(1), ulo(2):uhi(2), ulo(3):uhi(3),nvar)
  real(rt), intent(in) :: dx(3), prob_lo(3)
  
  integer :: i,j,k
  real(rt) :: x, Pt, rhot, uxt

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           x = prob_lo(1) + (i+0.5d0)*dx(1)

           if (x .lt. 0.5d0) then
              Pt = p_l
              rhot = rho_l
              uxt = u_l
           else
              Pt = p_r
              rhot = rho_r
              uxt = u_r
           end if

           u(i,j,k,urho) = rhot
           u(i,j,k,umx) = rhot*uxt
           u(i,j,k,umy:umz) = 0.d0
           u(i,j,k,ueint) = Pt / (gamma-1.d0)
           u(i,j,k,ueden) = u(i,j,k,ueint) + 0.5d0*rhot*uxt*uxt
           u(i,j,k,utemp) = u(i,j,k,ueint)/(u(i,j,k,urho)*cv)

        end do
     end do
  end do

end subroutine cns_initdata
