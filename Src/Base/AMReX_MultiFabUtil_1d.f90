! ***************************************************************************************
! subroutine bl_avgdown_with_vol
! ***************************************************************************************

subroutine bl_avgdown_with_vol (lo,hi,&
     fine,f_l1,f_h1, &
     crse,c_l1,c_h1, &
     fv,fv_l1,fv_h1, &
     lrat,ncomp)

  use amrex_fort_module, only : amrex_real
  implicit none

  integer f_l1,f_h1
  integer c_l1,c_h1
  integer fv_l1,fv_h1
  integer lo(1), hi(1)
  integer lrat(1), ncomp
  real(amrex_real) fine(f_l1:f_h1,ncomp)
  real(amrex_real) crse(c_l1:c_h1,ncomp)
  real(amrex_real) fv(fv_l1:fv_h1)

  integer :: i, ii, n, iref
  real(amrex_real) :: cv

  do n = 1, ncomp
     do i = lo(1), hi(1)
        ii = i * lrat(1)
        crse(i,n) = 0.d0
        cv        = 0.d0
        do iref = 0, lrat(1)-1
           cv        = cv        +                 fv(ii+iref)
           crse(i,n) = crse(i,n) + fine(ii+iref,n)*fv(ii+iref)
        end do
        crse(i,n) = crse(i,n) / cv
     end do
  end do

end subroutine bl_avgdown_with_vol


subroutine amrex_compute_divergence (lo, hi, divu, dlo, dhi, u, ulo, uhi, dxinv) bind(c)
  use amrex_fort_module, only : amrex_real
  implicit none
  integer, dimension(1), intent(in) :: lo, hi, dlo, dhi, ulo, uhi
  real(amrex_real), intent(inout) :: divu(dlo(1):dhi(1))
  real(amrex_real), intent(in   ) ::    u(ulo(1):uhi(1))
  real(amrex_real), intent(in) :: dxinv(1)
  integer :: i
  do i = lo(1), hi(1)
     divu(i) = dxinv(1) * (u(i+1)-u(i))
  end do
end subroutine amrex_compute_divergence
