
module amrex_mlnodelap_3d_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  use amrex_constants_module
  use amrex_lo_bctypes_module, only : amrex_lo_dirichlet, amrex_lo_neumann, amrex_lo_inflow, amrex_lo_periodic
  implicit none

  ! external dirichlet at physical boundary or internal dirichlet at crse/fine boundary
  integer, parameter :: dirichlet = 1

  integer, parameter :: crse_cell = 0
  integer, parameter :: fine_cell = 1
  integer, parameter :: crse_node = 0
  integer, parameter :: crse_fine_node = 1
  integer, parameter :: fine_node = 2

  private
  public :: &
       ! masks
       ! coeffs
       ! bc
       ! operator
       ! restriction
       ! interpolation
       ! rhs & u
       ! residual
       amrex_mlndlap_res_fine_contrib, amrex_mlndlap_res_cf_contrib
       ! sync residual

  ! RAP

contains

  subroutine amrex_mlndlap_res_fine_contrib (clo, chi, cglo, cghi, f, flo, fhi, &
       x, xlo, xhi, sig, slo, shi, Ax, alo, ahi, msk, mlo, mhi, dxinv) &
       bind(c,name='amrex_mlndlap_res_fine_contrib')
    integer, dimension(3), intent(in) :: clo, chi, cglo, cghi, flo, fhi, xlo, xhi, &
         slo, shi, alo, ahi, mlo, mhi
    real(amrex_real), intent(inout) :: f  (flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(amrex_real), intent(in   ) :: x  (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(amrex_real), intent(inout) :: Ax (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
    integer         , intent(in   ) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    real(amrex_real), intent(in) :: dxinv(3)

    integer, dimension(3) :: lo, hi, glo, ghi, gtlo, gthi
    integer :: i, j, k, ii, jj, kk, step
    real(amrex_real) :: facx, facy, facz, fxyz, fmx2y2z, f2xmy2z, f2x2ymz
    real(amrex_real) :: f4xm2ym2z, fm2x4ym2z, fm2xm2y4z
    real(amrex_real), parameter :: rfd = 0.125d0
    real(amrex_real), parameter :: chip = 0.5d0
    real(amrex_real), parameter :: chip2 = 0.25d0
    real(amrex_real), parameter :: chip3 = 0.125d0

    facx = (1.d0/36.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/36.d0)*dxinv(2)*dxinv(2)
    facz = (1.d0/36.d0)*dxinv(3)*dxinv(3)
    fxyz = facx + facy + facz
    fmx2y2z = -facx + 2.d0*facy + 2.d0*facz
    f2xmy2z = 2.d0*facx - facy + 2.d0*facz
    f2x2ymz = 2.d0*facx + 2.d0*facy - facz
    f4xm2ym2z = 4.d0*facx - 2.d0*facy - 2.d0*facz
    fm2x4ym2z = -2.d0*facx + 4.d0*facy - 2.d0*facz
    fm2xm2y4z = -2.d0*facx - 2.d0*facy + 4.d0*facz

    lo = 2*clo
    hi = 2*chi
    glo = 2*cglo
    ghi = 2*cghi

    gtlo = max(lo-1,glo)
    gthi = min(hi+1,ghi)

    do    kk = gtlo(3), gthi(3)
       do jj = gtlo(2), gthi(2)
          if (jj.eq.glo(2) .or. jj.eq.ghi(2) .or. kk.eq.glo(3) .or. kk.eq.ghi(3)) then
             step = 1
          else
             step = gthi(1)-gtlo(1)
          end if
          do ii = gtlo(1), gthi(1), step
             if (ii.eq.glo(1) .or. ii.eq.ghi(1) .or. step.eq.1) then
                Ax(ii,jj,kk) = x(ii,jj,kk)*(-4.d0)*fxyz* &
                     (sig(ii-1,jj-1,kk-1)+sig(ii,jj-1,kk-1)+sig(ii-1,jj,kk-1)+sig(ii,jj,kk-1) &
                     +sig(ii-1,jj-1,kk  )+sig(ii,jj-1,kk  )+sig(ii-1,jj,kk  )+sig(ii,jj,kk  )) &
                     !
                     + fxyz*(x(ii-1,jj-1,kk-1)*sig(ii-1,jj-1,kk-1) &
                     &     + x(ii+1,jj-1,kk-1)*sig(ii  ,jj-1,kk-1) &
                     &     + x(ii-1,jj+1,kk-1)*sig(ii-1,jj  ,kk-1) &
                     &     + x(ii+1,jj+1,kk-1)*sig(ii  ,jj  ,kk-1) &
                     &     + x(ii-1,jj-1,kk+1)*sig(ii-1,jj-1,kk  ) &
                     &     + x(ii+1,jj-1,kk+1)*sig(ii  ,jj-1,kk  ) &
                     &     + x(ii-1,jj+1,kk+1)*sig(ii-1,jj  ,kk  ) &
                     &     + x(ii+1,jj+1,kk+1)*sig(ii  ,jj  ,kk  )) &
                     !
                     + fmx2y2z*(x(ii  ,jj-1,kk-1)*(sig(ii-1,jj-1,kk-1)+sig(ii,jj-1,kk-1)) &
                     &        + x(ii  ,jj+1,kk-1)*(sig(ii-1,jj  ,kk-1)+sig(ii,jj  ,kk-1)) &
                     &        + x(ii  ,jj-1,kk+1)*(sig(ii-1,jj-1,kk  )+sig(ii,jj-1,kk  )) &
                     &        + x(ii  ,jj+1,kk+1)*(sig(ii-1,jj  ,kk  )+sig(ii,jj  ,kk  ))) &
                     !
                     + f2xmy2z*(x(ii-1,jj  ,kk-1)*(sig(ii-1,jj-1,kk-1)+sig(ii-1,jj,kk-1)) &
                     &        + x(ii+1,jj  ,kk-1)*(sig(ii  ,jj-1,kk-1)+sig(ii  ,jj,kk-1)) &
                     &        + x(ii-1,jj  ,kk+1)*(sig(ii-1,jj-1,kk  )+sig(ii-1,jj,kk  )) &
                     &        + x(ii+1,jj  ,kk+1)*(sig(ii  ,jj-1,kk  )+sig(ii  ,jj,kk  ))) &
                     !
                     + f2x2ymz*(x(ii-1,jj-1,kk  )*(sig(ii-1,jj-1,kk-1)+sig(ii-1,jj-1,kk)) &
                     &        + x(ii+1,jj-1,kk  )*(sig(ii  ,jj-1,kk-1)+sig(ii  ,jj-1,kk)) &
                     &        + x(ii-1,jj+1,kk  )*(sig(ii-1,jj  ,kk-1)+sig(ii-1,jj  ,kk)) &
                     &        + x(ii+1,jj+1,kk  )*(sig(ii  ,jj  ,kk-1)+sig(ii  ,jj  ,kk))) &
                     !
                     + f4xm2ym2z*(x(ii-1,jj,kk)*(sig(ii-1,jj-1,kk-1)+sig(ii-1,jj,kk-1)+sig(ii-1,jj-1,kk)+sig(ii-1,jj,kk)) &
                     &          + x(ii+1,jj,kk)*(sig(ii  ,jj-1,kk-1)+sig(ii  ,jj,kk-1)+sig(ii  ,jj-1,kk)+sig(ii  ,jj,kk))) &
                     + fm2x4ym2z*(x(ii,jj-1,kk)*(sig(ii-1,jj-1,kk-1)+sig(ii,jj-1,kk-1)+sig(ii-1,jj-1,kk)+sig(ii,jj-1,kk)) &
                     &          + x(ii,jj+1,kk)*(sig(ii-1,jj  ,kk-1)+sig(ii,jj  ,kk-1)+sig(ii-1,jj  ,kk)+sig(ii,jj  ,kk))) &
                     + fm2xm2y4z*(x(ii,jj,kk-1)*(sig(ii-1,jj-1,kk-1)+sig(ii,jj-1,kk-1)+sig(ii-1,jj,kk-1)+sig(ii,jj,kk-1)) &
                     &          + x(ii,jj,kk+1)*(sig(ii-1,jj-1,kk  )+sig(ii,jj-1,kk  )+sig(ii-1,jj,kk  )+sig(ii,jj,kk  )))

             end if
          end do
       end do
    end do

    do k = clo(3), chi(3)
       kk = 2*k
       do j = clo(2), chi(2)
          jj = 2*j
          if (jj.eq.glo(2) .or. jj.eq.ghi(2) .or. kk.eq.glo(3) .or. kk.eq.ghi(3)) then
             step = 1
          else
             step = chi(1)-clo(1)
          end if
          do i = clo(1), chi(1), step
             ii = 2*i
             if (msk(ii,jj,kk) .eq. dirichlet) then
                f(i,j,k) = f(i,j,k) + rfd*(Ax(ii,jj,kk) &
                     + chip*(Ax(ii,jj,kk-1)+Ax(ii,jj,kk+1) &
                     &      +Ax(ii,jj-1,kk)+Ax(ii,jj+1,kk) &
                     &      +Ax(ii-1,jj,kk)+Ax(ii+1,jj,kk)) &
                     + chip2*(Ax(ii,jj-1,kk-1)+Ax(ii,jj+1,kk-1)+Ax(ii,jj-1,kk+1)+Ax(ii,jj+1,kk+1) &
                     &       +Ax(ii-1,jj,kk-1)+Ax(ii+1,jj,kk-1)+Ax(ii-1,jj,kk+1)+Ax(ii+1,jj,kk+1) &
                     &       +Ax(ii-1,jj-1,kk)+Ax(ii+1,jj-1,kk)+Ax(ii-1,jj+1,kk)+Ax(ii+1,jj+1,kk)) &
                     + chip3*(Ax(ii-1,jj-1,kk-1)+Ax(ii+1,jj-1,kk-1) &
                     &       +Ax(ii-1,jj+1,kk-1)+Ax(ii+1,jj+1,kk-1) &
                     &       +Ax(ii-1,jj-1,kk+1)+Ax(ii+1,jj-1,kk+1) &
                     &       +Ax(ii-1,jj+1,kk+1)+Ax(ii+1,jj+1,kk+1)))
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_res_fine_contrib


  subroutine amrex_mlndlap_res_cf_contrib (lo, hi, res, rlo, rhi, phi, phlo, phhi, &
       rhs, rhlo, rhhi, sig, slo, shi, dmsk, mlo, mhi, ndmsk, nmlo, nmhi, ccmsk, cmlo, cmhi, &
       fc, clo, chi, dxinv, ndlo, ndhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_res_cf_contrib')
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, phlo, phhi, rhlo, rhhi, slo, shi, &
         mlo, mhi, nmlo, nmhi, cmlo, cmhi, clo, chi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: res( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) :: phi(phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))
    real(amrex_real), intent(in   ) :: rhs(rhlo(1):rhhi(1),rhlo(2):rhhi(2),rhlo(3):rhhi(3))
    real(amrex_real), intent(in   ) :: sig( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3))
    real(amrex_real), intent(inout) :: fc ( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3))
    integer, intent(in) ::  dmsk( mlo(1): mhi(1), mlo(2): mhi(2), mlo(3): mhi(3))
    integer, intent(in) :: ndmsk(nmlo(1):nmhi(1),nmlo(2):nmhi(2),nmlo(3):nmhi(3))
    integer, intent(in) :: ccmsk(cmlo(1):cmhi(1),cmlo(2):cmhi(2),cmlo(3):cmhi(3))

    integer :: i,j,k
    real(amrex_real) :: Ax, Axf, facx, facy, facz

    facx = (1.d0/36.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/36.d0)*dxinv(2)*dxinv(2)
    facz = (1.d0/36.d0)*dxinv(3)*dxinv(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (dmsk(i,j,k) .ne. dirichlet) then
                if (ndmsk(i,j,k) .eq. crse_fine_node) then
                   Ax = 0.d0
                   if (ccmsk(i-1,j-1,k-1) .eq. crse_cell) then
                      Ax = Ax + sig(i-1,j-1,k-1)*(facx*(4.d0*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  )) &
                           &                           +2.d0*(phi(i-1,j-1,k  )-phi(i  ,j-1,k  )) &
                           &                           +2.d0*(phi(i-1,j  ,k-1)-phi(i  ,j  ,k-1)) &
                           &                           +     (phi(i-1,j-1,k-1)-phi(i  ,j-1,k-1))) &
                           &                    + facy*(4.d0*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  )) &
                           &                           +2.d0*(phi(i-1,j-1,k  )-phi(i-1,j  ,k  )) &
                           &                           +2.d0*(phi(i  ,j-1,k-1)-phi(i  ,j  ,k-1)) &
                           &                           +     (phi(i-1,j-1,k-1)-phi(i-1,j  ,k-1))) &
                           &                    + facz*(4.d0*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  )) &
                           &                           +2.d0*(phi(i-1,j  ,k-1)-phi(i-1,j  ,k  )) &
                           &                           +2.d0*(phi(i  ,j-1,k-1)-phi(i  ,j-1,k  )) &
                           &                           +     (phi(i-1,j-1,k-1)-phi(i-1,j-1,k  ))))
                   end if
                   if (ccmsk(i,j-1,k-1) .eq. crse_cell) then
                      Ax = Ax + sig(i,j-1,k-1)*(facx*(4.d0*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i+1,j-1,k  )-phi(i  ,j-1,k  )) &
                           &                         +2.d0*(phi(i+1,j  ,k-1)-phi(i  ,j  ,k-1)) &
                           &                         +     (phi(i+1,j-1,k-1)-phi(i  ,j-1,k-1))) &
                           &                  + facy*(4.d0*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i+1,j-1,k  )-phi(i+1,j  ,k  )) &
                           &                         +2.d0*(phi(i  ,j-1,k-1)-phi(i  ,j  ,k-1)) &
                           &                         +     (phi(i+1,j-1,k-1)-phi(i+1,j  ,k-1))) &
                           &                  + facz*(4.d0*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i+1,j  ,k-1)-phi(i+1,j  ,k  )) &
                           &                         +2.d0*(phi(i  ,j-1,k-1)-phi(i  ,j-1,k  )) &
                           &                         +     (phi(i+1,j-1,k-1)-phi(i+1,j-1,k  ))))
                   end if
                   if (ccmsk(i-1,j,k-1) .eq. crse_cell) then
                      Ax = Ax + sig(i-1,j,k-1)*(facx*(4.d0*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i-1,j+1,k  )-phi(i  ,j+1,k  )) &
                           &                         +2.d0*(phi(i-1,j  ,k-1)-phi(i  ,j  ,k-1)) &
                           &                         +     (phi(i-1,j+1,k-1)-phi(i  ,j+1,k-1))) &
                           &                  + facy*(4.d0*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i-1,j+1,k  )-phi(i-1,j  ,k  )) &
                           &                         +2.d0*(phi(i  ,j+1,k-1)-phi(i  ,j  ,k-1)) &
                           &                         +     (phi(i-1,j+1,k-1)-phi(i-1,j  ,k-1))) &
                           &                  + facz*(4.d0*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i-1,j  ,k-1)-phi(i-1,j  ,k  )) &
                           &                         +2.d0*(phi(i  ,j+1,k-1)-phi(i  ,j+1,k  )) &
                           &                         +     (phi(i-1,j+1,k-1)-phi(i-1,j+1,k  ))))
                   end if
                   if (ccmsk(i,j,k-1) .eq. crse_cell) then
                      Ax = Ax + sig(i,j,k-1)*(facx*(4.d0*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i+1,j+1,k  )-phi(i  ,j+1,k  )) &
                           &                       +2.d0*(phi(i+1,j  ,k-1)-phi(i  ,j  ,k-1)) &
                           &                       +     (phi(i+1,j+1,k-1)-phi(i  ,j+1,k-1))) &
                           &                + facy*(4.d0*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i+1,j+1,k  )-phi(i+1,j  ,k  )) &
                           &                       +2.d0*(phi(i  ,j+1,k-1)-phi(i  ,j  ,k-1)) &
                           &                       +     (phi(i+1,j+1,k-1)-phi(i+1,j  ,k-1))) &
                           &                + facz*(4.d0*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i+1,j  ,k-1)-phi(i+1,j  ,k  )) &
                           &                       +2.d0*(phi(i  ,j+1,k-1)-phi(i  ,j+1,k  )) &
                           &                       +     (phi(i+1,j+1,k-1)-phi(i+1,j+1,k  ))))
                   end if
                   if (ccmsk(i-1,j-1,k) .eq. crse_cell) then
                      Ax = Ax + sig(i-1,j-1,k)*(facx*(4.d0*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i-1,j-1,k  )-phi(i  ,j-1,k  )) &
                           &                         +2.d0*(phi(i-1,j  ,k+1)-phi(i  ,j  ,k+1)) &
                           &                         +     (phi(i-1,j-1,k+1)-phi(i  ,j-1,k+1))) &
                           &                  + facy*(4.d0*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i-1,j-1,k  )-phi(i-1,j  ,k  )) &
                           &                         +2.d0*(phi(i  ,j-1,k+1)-phi(i  ,j  ,k+1)) &
                           &                         +     (phi(i-1,j-1,k+1)-phi(i-1,j  ,k+1))) &
                           &                  + facz*(4.d0*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i-1,j  ,k+1)-phi(i-1,j  ,k  )) &
                           &                         +2.d0*(phi(i  ,j-1,k+1)-phi(i  ,j-1,k  )) &
                           &                         +     (phi(i-1,j-1,k+1)-phi(i-1,j-1,k  ))))
                   end if
                   if (ccmsk(i,j-1,k) .eq. crse_cell) then
                      Ax = Ax + sig(i,j-1,k)*(facx*(4.d0*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i+1,j-1,k  )-phi(i  ,j-1,k  )) &
                           &                       +2.d0*(phi(i+1,j  ,k+1)-phi(i  ,j  ,k+1)) &
                           &                       +     (phi(i+1,j-1,k+1)-phi(i  ,j-1,k+1))) &
                           &                + facy*(4.d0*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i+1,j-1,k  )-phi(i+1,j  ,k  )) &
                           &                       +2.d0*(phi(i  ,j-1,k+1)-phi(i  ,j  ,k+1)) &
                           &                       +     (phi(i+1,j-1,k+1)-phi(i+1,j  ,k+1))) &
                           &                + facz*(4.d0*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i+1,j  ,k+1)-phi(i+1,j  ,k  )) &
                           &                       +2.d0*(phi(i  ,j-1,k+1)-phi(i  ,j-1,k  )) &
                           &                       +     (phi(i+1,j-1,k+1)-phi(i+1,j-1,k  ))))
                   end if
                   if (ccmsk(i-1,j,k) .eq. crse_cell) then
                      Ax = Ax + sig(i-1,j,k)*(facx*(4.d0*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i-1,j+1,k  )-phi(i  ,j+1,k  )) &
                           &                       +2.d0*(phi(i-1,j  ,k+1)-phi(i  ,j  ,k+1)) &
                           &                       +     (phi(i-1,j+1,k+1)-phi(i  ,j+1,k+1))) &
                           &                + facy*(4.d0*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i-1,j+1,k  )-phi(i-1,j  ,k  )) &
                           &                       +2.d0*(phi(i  ,j+1,k+1)-phi(i  ,j  ,k+1)) &
                           &                       +     (phi(i-1,j+1,k+1)-phi(i-1,j  ,k+1))) &
                           &                + facz*(4.d0*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i-1,j  ,k+1)-phi(i-1,j  ,k  )) &
                           &                       +2.d0*(phi(i  ,j+1,k+1)-phi(i  ,j+1,k  )) &
                           &                       +     (phi(i-1,j+1,k+1)-phi(i-1,j+1,k  ))))
                   end if
                   if (ccmsk(i,j,k) .eq. crse_cell) then
                      Ax = Ax + sig(i,j,k)*(facx*(4.d0*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  )) &
                           &                     +2.d0*(phi(i+1,j+1,k  )-phi(i  ,j+1,k  )) &
                           &                     +2.d0*(phi(i+1,j  ,k+1)-phi(i  ,j  ,k+1)) &
                           &                     +     (phi(i+1,j+1,k+1)-phi(i  ,j+1,k+1))) &
                           &              + facy*(4.d0*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  )) &
                           &                     +2.d0*(phi(i+1,j+1,k  )-phi(i+1,j  ,k  )) &
                           &                     +2.d0*(phi(i  ,j+1,k+1)-phi(i  ,j  ,k+1)) &
                           &                     +     (phi(i+1,j+1,k+1)-phi(i+1,j  ,k+1))) &
                           &              + facz*(4.d0*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  )) &
                           &                     +2.d0*(phi(i+1,j  ,k+1)-phi(i+1,j  ,k  )) &
                           &                     +2.d0*(phi(i  ,j+1,k+1)-phi(i  ,j+1,k  )) &
                           &                     +     (phi(i+1,j+1,k+1)-phi(i+1,j+1,k  ))))
                   end if

                   Axf = fc(i,j,k)

                   if (i .eq. ndlo(1) .and. &
                        (    bclo(1) .eq. amrex_lo_neumann &
                        .or. bclo(1) .eq. amrex_lo_inflow)) then
                      Axf = 2.d0*Axf
                   else if (i.eq. ndhi(1) .and. &
                        (    bchi(1) .eq. amrex_lo_neumann &
                        .or. bchi(1) .eq. amrex_lo_inflow)) then
                      Axf = 2.d0*Axf
                   end if

                   if (j .eq. ndlo(2) .and. &
                        (    bclo(2) .eq. amrex_lo_neumann &
                        .or. bclo(2) .eq. amrex_lo_inflow)) then
                      Axf = 2.d0*Axf
                   else if (j .eq. ndhi(2) .and. &
                        (    bchi(2) .eq. amrex_lo_neumann &
                        .or. bchi(2) .eq. amrex_lo_inflow)) then
                      Axf = 2.d0*Axf
                   end if

                   if (k .eq. ndlo(3) .and. &
                        (    bclo(3) .eq. amrex_lo_neumann &
                        .or. bclo(3) .eq. amrex_lo_inflow)) then
                      Axf = 2.d0*Axf
                   else if (k .eq. ndhi(3) .and. &
                        (    bchi(3) .eq. amrex_lo_neumann &
                        .or. bchi(3) .eq. amrex_lo_inflow)) then
                      Axf = 2.d0*Axf
                   end if

                   res(i,j,k) = rhs(i,j,k) - (Ax + Axf)
                end if
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_res_cf_contrib


end module amrex_mlnodelap_3d_module
