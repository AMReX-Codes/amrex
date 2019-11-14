module amrex_mlnodelap_2d_module

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

  logical, private, save :: is_rz = .false.

  private
  public :: &
       amrex_mlndlap_set_rz, &
       ! masks
       amrex_mlndlap_any_fine_sync_cells, &
       ! coeffs
       ! bc
       ! operator
       ! restriction
       ! interpolation
       ! rhs & u
       ! residual
       amrex_mlndlap_crse_resid, &
       amrex_mlndlap_res_fine_contrib, amrex_mlndlap_res_cf_contrib, &
       ! sync residual
       amrex_mlndlap_zero_fine

contains

  subroutine amrex_mlndlap_set_rz (rz) bind(c,name='amrex_mlndlap_set_rz')
    integer, intent(in) :: rz
    is_rz = rz.ne.0
  end subroutine amrex_mlndlap_set_rz


  function amrex_mlndlap_any_fine_sync_cells (lo, hi, msk, mlo, mhi, fine_flag) result(r) &
       bind(c,name='amrex_mlndlap_any_fine_sync_cells')
    integer :: r
    integer, dimension(2), intent(in) :: lo, hi, mlo, mhi
    integer, intent(in   ) :: msk  ( mlo(1): mhi(1), mlo(2): mhi(2))
    integer, intent(in) :: fine_flag

    integer :: i,j

    r = 0

    do j = lo(2), hi(2)
       if (r.eq.1) exit
       do i = lo(1), hi(1)
          if (r.eq.1) exit
          if (msk(i,j) .eq. fine_flag) r = 1
       end do
    end do
  end function amrex_mlndlap_any_fine_sync_cells


  subroutine amrex_mlndlap_crse_resid (lo, hi, resid, rslo, rshi, rhs, rhlo, rhhi, msk, mlo, mhi, &
       ndlo, ndhi, bclo, bchi) bind(c, name='amrex_mlndlap_crse_resid')
    integer, dimension(2), intent(in) :: lo, hi, rslo, rshi, rhlo, rhhi, mlo, mhi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(inout) :: resid(rslo(1):rshi(1),rslo(2):rshi(2))
    real(amrex_real), intent(in   ) :: rhs  (rhlo(1):rhhi(1),rhlo(2):rhhi(2))
    integer         , intent(in   ) :: msk  ( mlo(1): mhi(1), mlo(2): mhi(2))

    integer :: i,j
    real(amrex_real) :: fac

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (any(msk(i-1:i,j-1:j).eq.0) .and. any(msk(i-1:i,j-1:j).eq.1)) then

             fac = 1.d0
             if (i .eq. ndlo(1) .and. &
                  (    bclo(1) .eq. amrex_lo_neumann &
                  .or. bclo(1) .eq. amrex_lo_inflow)) then
                fac = fac*2.d0
             else if (i.eq. ndhi(1) .and. &
                  (    bchi(1) .eq. amrex_lo_neumann &
                  .or. bchi(1) .eq. amrex_lo_inflow)) then
                fac = fac*2.d0
             end if
             
             if (j .eq. ndlo(2) .and. &
                  (    bclo(2) .eq. amrex_lo_neumann &
                  .or. bclo(2) .eq. amrex_lo_inflow)) then
                fac = fac*2.d0
             else if (j .eq. ndhi(2) .and. &
                  (    bchi(2) .eq. amrex_lo_neumann &
                  .or. bchi(2) .eq. amrex_lo_inflow)) then
                fac = fac*2.d0
             end if

             resid(i,j) = (rhs(i,j) - resid(i,j)) * fac
          else
             resid(i,j) = 0.d0
          end if
       end do
    end do
  end subroutine amrex_mlndlap_crse_resid


  subroutine amrex_mlndlap_res_fine_contrib (clo, chi, cglo, cghi, f, flo, fhi, &
       x, xlo, xhi, sig, slo, shi, Ax, alo, ahi, msk, mlo, mhi, dxinv) &
       bind(c,name='amrex_mlndlap_res_fine_contrib')
    integer, dimension(2), intent(in) :: clo, chi, cglo, cghi, flo, fhi, xlo, xhi, &
         slo, shi, alo, ahi, mlo, mhi
    real(amrex_real), intent(inout) :: f  (flo(1):fhi(1),flo(2):fhi(2))
    real(amrex_real), intent(in   ) :: x  (xlo(1):xhi(1),xlo(2):xhi(2))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1),slo(2):shi(2))
    real(amrex_real), intent(inout) :: Ax (alo(1):ahi(1),alo(2):ahi(2))
    integer         , intent(in   ) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))
    real(amrex_real), intent(in) :: dxinv(2)

    integer, dimension(2) :: lo, hi, glo, ghi, gtlo, gthi
    integer :: i, j, ii, jj, step
    real(amrex_real) :: facx, facy, fxy, f2xmy, fmx2y, fm, fp
    real(amrex_real), parameter :: rfd = 0.25d0
    real(amrex_real), parameter :: chip = 0.5d0
    real(amrex_real), parameter :: chip2 = 0.25d0

    facx = (1.d0/6.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/6.d0)*dxinv(2)*dxinv(2)
    fxy = facx + facy
    f2xmy = 2.d0*facx - facy
    fmx2y = 2.d0*facy - facx

    lo = 2*clo
    hi = 2*chi
    glo = 2*cglo
    ghi = 2*cghi

    gtlo = max(lo-1,glo)
    gthi = min(hi+1,ghi)

    do jj = gtlo(2), gthi(2)
       if (jj .eq. glo(2) .or. jj .eq. ghi(2)) then
          step = 1
       else
          step = gthi(1)-gtlo(1)
       end if
       do ii = gtlo(1), gthi(1), step
          if (ii.eq.glo(1) .or. ii.eq.ghi(1) .or. step .eq. 1) then
             Ax(ii,jj) = x(ii-1,jj-1)*fxy*sig(ii-1,jj-1) &
                  +      x(ii+1,jj-1)*fxy*sig(ii  ,jj-1) &
                  +      x(ii-1,jj+1)*fxy*sig(ii-1,jj  ) &
                  +      x(ii+1,jj+1)*fxy*sig(ii  ,jj  ) &
                  +      x(ii-1,jj)*f2xmy*(sig(ii-1,jj-1)+sig(ii-1,jj  )) &
                  +      x(ii+1,jj)*f2xmy*(sig(ii  ,jj-1)+sig(ii  ,jj  )) &
                  +      x(ii,jj-1)*fmx2y*(sig(ii-1,jj-1)+sig(ii  ,jj-1)) &
                  +      x(ii,jj+1)*fmx2y*(sig(ii-1,jj  )+sig(ii  ,jj  )) &
                  +      x(ii,jj)*(-2.d0)*fxy*(sig(ii-1,jj-1)+sig(ii,jj-1)+sig(ii-1,jj)+sig(ii,jj))

             if (is_rz) then
                fp = facy / (2*ii+1)
                fm = facy / (2*ii-1)
                Ax(ii,jj) = Ax(ii,jj) + (fm*sig(ii-1,jj  )-fp*sig(ii,jj  ))*(x(ii,jj+1)-x(ii,jj)) &
                     &                + (fm*sig(ii-1,jj-1)-fp*sig(ii,jj-1))*(x(ii,jj-1)-x(ii,jj))
             end if
          end if
       end do
    end do

    do j = clo(2), chi(2)
       jj = 2*j
       if (j .eq. cglo(2) .or. j .eq. cghi(2)) then
          step = 1
       else
          step = chi(1)-clo(1)
       end if
       do i = clo(1), chi(1), step
          ii = 2*i
          if (msk(ii,jj) .eq. dirichlet) then
             f(i,j) = f(i,j) + rfd*(Ax(ii,jj) &
                  + chip*(Ax(ii-1,jj)+Ax(ii+1,jj)+Ax(ii,jj-1)+Ax(ii,jj+1)) &
                  + chip2*(Ax(ii-1,jj-1)+Ax(ii+1,jj-1)+Ax(ii-1,jj+1)+Ax(ii+1,jj+1)))
          end if
       end do
    end do

  end subroutine amrex_mlndlap_res_fine_contrib


  subroutine amrex_mlndlap_res_cf_contrib (lo, hi, res, rlo, rhi, phi, phlo, phhi, &
       rhs, rhlo, rhhi, sig, slo, shi, dmsk, mlo, mhi, ndmsk, nmlo, nmhi, ccmsk, cmlo, cmhi, &
       fc, clo, chi, dxinv, ndlo, ndhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_res_cf_contrib')
    integer, dimension(2), intent(in) :: lo, hi, rlo, rhi, phlo, phhi, rhlo, rhhi, slo, shi, &
         mlo, mhi, nmlo, nmhi, cmlo, cmhi, clo, chi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: res( rlo(1): rhi(1), rlo(2): rhi(2))
    real(amrex_real), intent(in   ) :: phi(phlo(1):phhi(1),phlo(2):phhi(2))
    real(amrex_real), intent(in   ) :: rhs(rhlo(1):rhhi(1),rhlo(2):rhhi(2))
    real(amrex_real), intent(in   ) :: sig( slo(1): shi(1), slo(2): shi(2))
    real(amrex_real), intent(inout) :: fc ( clo(1): chi(1), clo(2): chi(2))
    integer, intent(in) ::  dmsk( mlo(1): mhi(1), mlo(2): mhi(2))
    integer, intent(in) :: ndmsk(nmlo(1):nmhi(1),nmlo(2):nmhi(2))
    integer, intent(in) :: ccmsk(cmlo(1):cmhi(1),cmlo(2):cmhi(2))

    integer :: i,j
    real(amrex_real) :: Ax, Axf, facx, facy, fp, fm

    facx = (1.d0/6.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/6.d0)*dxinv(2)*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (dmsk(i,j) .ne. dirichlet) then
             if (ndmsk(i,j) .eq. crse_fine_node) then

                if (is_rz) then
                   fp = facy / (2*i+1)
                   fm = facy / (2*i-1)
                end if

                Ax = 0.d0
                if (ccmsk(i-1,j-1) .eq. crse_cell) then
                   Ax = Ax + sig(i-1,j-1)*(facx*(2.d0*(phi(i-1,j  )-phi(i  ,j  )) &
                        &                       +     (phi(i-1,j-1)-phi(i  ,j-1))) &
                        &                + facy*(2.d0*(phi(i  ,j-1)-phi(i  ,j  )) &
                        &                       +     (phi(i-1,j-1)-phi(i-1,j  ))))
                   if (is_rz) then
                      Ax = Ax + fm*sig(i-1,j-1)*(phi(i,j-1)-phi(i,j))
                   end if
                end if
                if (ccmsk(i,j-1) .eq. crse_cell) then
                   Ax = Ax + sig(i,j-1)*(facx*(2.d0*(phi(i+1,j  )-phi(i  ,j  )) &
                        &                     +     (phi(i+1,j-1)-phi(i  ,j-1))) &
                        &              + facy*(2.d0*(phi(i  ,j-1)-phi(i  ,j  )) &
                        &                     +     (phi(i+1,j-1)-phi(i+1,j  ))))
                   if (is_rz) then
                      Ax = Ax - fp*sig(i,j-1)*(phi(i,j-1)-phi(i,j))
                   end if
                end if
                if (ccmsk(i-1,j) .eq. crse_cell) then
                   Ax = Ax + sig(i-1,j)*(facx*(2.d0*(phi(i-1,j  )-phi(i  ,j  )) &
                        &                     +     (phi(i-1,j+1)-phi(i  ,j+1))) &
                        &              + facy*(2.d0*(phi(i  ,j+1)-phi(i  ,j  )) &
                        &                     +     (phi(i-1,j+1)-phi(i-1,j  ))))
                   if (is_rz) then
                      Ax = Ax + fm*sig(i-1,j)*(phi(i,j+1)-phi(i,j))
                   end if
                end if
                if (ccmsk(i,j) .eq. crse_cell) then
                   Ax = Ax + sig(i,j)*(facx*(2.d0*(phi(i+1,j  )-phi(i  ,j  )) &
                        &                  +      (phi(i+1,j+1)-phi(i  ,j+1))) &
                        &            + facy*(2.d0*(phi(i  ,j+1)-phi(i  ,j  )) &
                        &                  +      (phi(i+1,j+1)-phi(i+1,j  ))))
                   if (is_rz) then
                      Ax = Ax - fp*sig(i,j)*(phi(i,j+1)-phi(i,j))
                   end if
                end if

                Axf = fc(i,j)

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

                res(i,j) = rhs(i,j) - (Ax + Axf)
             end if
          end if
       end do
    end do
  end subroutine amrex_mlndlap_res_cf_contrib


  subroutine amrex_mlndlap_zero_fine (lo, hi, phi, dlo, dhi, msk, mlo, mhi, fine_flag) &
       bind(c, name='amrex_mlndlap_zero_fine')
    integer, dimension(2), intent(in) :: lo, hi, dlo, dhi, mlo, mhi
    real(amrex_real), intent(inout) :: phi(dlo(1):dhi(1),dlo(2):dhi(2))
    integer         , intent(in   ) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))
    integer, intent(in) :: fine_flag

    integer :: i,j

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          ! Testing if the node is covered by a fine level in computing
          ! coarse sync residual
          if (all(msk(i-1:i,j-1:j).eq.fine_flag)) then
             phi(i,j) = 0.d0
          end if
       end do
    end do
  end subroutine amrex_mlndlap_zero_fine

end module amrex_mlnodelap_2d_module
