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

  real(amrex_real), private, parameter :: eps = 1.d-100

#ifdef AMREX_USE_EB
  integer, private, parameter :: i_S_x     = 1
  integer, private, parameter :: i_S_y     = 2
  integer, private, parameter :: i_S_x2    = 3
  integer, private, parameter :: i_S_y2    = 4
  integer, private, parameter :: i_S_xy    = 5
  integer, private, parameter :: n_Sintg   = 5
#endif

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
       amrex_mlndlap_divu_fine_contrib, amrex_mlndlap_divu_cf_contrib, &
       amrex_mlndlap_rhcc_fine_contrib, amrex_mlndlap_rhcc_crse_contrib, &
       ! residual
       amrex_mlndlap_crse_resid, &
       amrex_mlndlap_res_fine_contrib, amrex_mlndlap_res_cf_contrib, &
       ! sync residual
       amrex_mlndlap_zero_fine

  ! RAP
  public:: &
       amrex_mlndlap_stencil_rap

#ifdef AMREX_USE_EB
  public:: amrex_mlndlap_set_integral, amrex_mlndlap_set_integral_eb, &
       amrex_mlndlap_set_connection, amrex_mlndlap_set_stencil_eb, &
       amrex_mlndlap_divu_eb, amrex_mlndlap_mknewu_eb, amrex_mlndlap_rhcc_eb
#endif

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


  subroutine amrex_mlndlap_divu_fine_contrib (clo, chi, cglo, cghi, rhs, rlo, rhi, &
       vel, vlo, vhi, frh, flo, fhi, msk, mlo, mhi, dxinv) &
       bind(c,name='amrex_mlndlap_divu_fine_contrib')
    integer, dimension(2), intent(in) :: clo, chi, cglo, cghi, rlo, rhi, vlo, vhi, &
         flo, fhi, mlo, mhi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1),vlo(2):vhi(2),2)
    real(amrex_real), intent(inout) :: frh(flo(1):fhi(1),flo(2):fhi(2))
    integer         , intent(in   ) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))

    integer, dimension(2) :: lo, hi, glo, ghi, gtlo, gthi
    integer :: i, j, ii, jj, step
    real(amrex_real) :: facx, facy, fm, fp
    real(amrex_real), parameter :: rfd = 0.25d0
    real(amrex_real), parameter :: chip = 0.5d0
    real(amrex_real), parameter :: chip2 = 0.25d0

    ! note that dxinv is fine dxinv
    facx = 0.5d0*dxinv(1)
    facy = 0.5d0*dxinv(2)

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
             frh(ii,jj) = facx*(-vel(ii-1,jj-1,1)+vel(ii,jj-1,1)-vel(ii-1,jj,1)+vel(ii,jj,1)) &
                    &   + facy*(-vel(ii-1,jj-1,2)-vel(ii,jj-1,2)+vel(ii-1,jj,2)+vel(ii,jj,2))

             if (is_rz) then
                fm = facy / (6*ii-3)
                fp = facy / (6*ii+3)
                frh(ii,jj) = frh(ii,jj) + fm*(vel(ii-1,jj,2)-vel(ii-1,jj-1,2)) &
                     &                  - fp*(vel(ii  ,jj,2)-vel(ii  ,jj-1,2))
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
             rhs(i,j) = rhs(i,j) + rfd*(frh(ii,jj) &
                  + chip*(frh(ii-1,jj)+frh(ii+1,jj)+frh(ii,jj-1)+frh(ii,jj+1)) &
                  + chip2*(frh(ii-1,jj-1)+frh(ii+1,jj-1)+frh(ii-1,jj+1)+frh(ii+1,jj+1)))
          end if
       end do
    end do

  end subroutine amrex_mlndlap_divu_fine_contrib


  subroutine amrex_mlndlap_divu_cf_contrib (lo, hi,  rhs, rlo, rhi, vel, vlo, vhi, dmsk, mlo, mhi, &
       ndmsk, nmlo, nmhi, ccmsk, cmlo, cmhi, fc, clo, chi, dxinv, ndlo, ndhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_divu_cf_contrib')
    integer, dimension(2), intent(in) :: lo, hi, rlo, rhi, vlo, vhi, mlo, mhi, &
         nmlo, nmhi, cmlo, cmhi, clo, chi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1),vlo(2):vhi(2),2)
    real(amrex_real), intent(in   ) :: fc (clo(1):chi(1),clo(2):chi(2))
    integer, intent(in) :: dmsk(mlo(1):mhi(1),mlo(2):mhi(2))
    integer, intent(in) :: ndmsk(nmlo(1):nmhi(1),nmlo(2):nmhi(2))
    integer, intent(in) :: ccmsk(cmlo(1):cmhi(1),cmlo(2):cmhi(2))

    integer :: i,j
    real(amrex_real) :: facx, facy, fm, fp

    facx = 0.5d0*dxinv(1)
    facy = 0.5d0*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (dmsk(i,j) .ne. dirichlet) then
             if (ndmsk(i,j) .eq. crse_fine_node) then
                rhs(i,j) = fc(i,j) &
                     + (1.d0-ccmsk(i-1,j-1)) * (-facx*vel(i-1,j-1,1) - facy*vel(i-1,j-1,2)) &
                     + (1.d0-ccmsk(i  ,j-1)) * ( facx*vel(i  ,j-1,1) - facy*vel(i  ,j-1,2)) &
                     + (1.d0-ccmsk(i-1,j  )) * (-facx*vel(i-1,j  ,1) + facy*vel(i-1,j  ,2)) &
                     + (1.d0-ccmsk(i  ,j  )) * ( facx*vel(i  ,j  ,1) + facy*vel(i  ,j  ,2))

                if (is_rz) then
                   fm = facy / (6*i-3)
                   fp = facy / (6*i+3)
                   rhs(i,j) = rhs(i,j) + fm*((1.d0-ccmsk(i-1,j  ))*vel(i-1,j  ,2) &
                        &                   -(1.d0-ccmsk(i-1,j-1))*vel(i-1,j-1,2)) &
                        &              - fp*((1.d0-ccmsk(i  ,j  ))*vel(i  ,j  ,2) &
                        &                   -(1.d0-ccmsk(i  ,j-1))*vel(i  ,j-1,2))
                end if

                if (i .eq. ndlo(1) .and. &
                     (    bclo(1) .eq. amrex_lo_neumann &
                     .or. bclo(1) .eq. amrex_lo_inflow)) then
                   rhs(i,j) = 2.d0*rhs(i,j)
                else if (i.eq. ndhi(1) .and. &
                     (    bchi(1) .eq. amrex_lo_neumann &
                     .or. bchi(1) .eq. amrex_lo_inflow)) then
                   rhs(i,j) = 2.d0*rhs(i,j)
                end if

                if (j .eq. ndlo(2) .and. &
                     (    bclo(2) .eq. amrex_lo_neumann &
                     .or. bclo(2) .eq. amrex_lo_inflow)) then
                   rhs(i,j) = 2.d0*rhs(i,j)                   
                else if (j .eq. ndhi(2) .and. &
                     (    bchi(2) .eq. amrex_lo_neumann &
                     .or. bchi(2) .eq. amrex_lo_inflow)) then
                   rhs(i,j) = 2.d0*rhs(i,j)
                end if
             end if
          end if
       end do
    end do

  end subroutine amrex_mlndlap_divu_cf_contrib


  subroutine amrex_mlndlap_rhcc_fine_contrib (clo, chi, cglo, cghi, rhs, rlo, rhi, &
       cc, cclo, cchi, msk, mlo, mhi) bind(c,name='amrex_mlndlap_rhcc_fine_contrib')
    integer, dimension(2), intent(in) :: clo, chi, cglo, cghi, rlo, rhi, cclo, cchi, mlo, mhi
    real(amrex_real), intent(inout) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2))
    real(amrex_real), intent(in   ) :: cc (cclo(1):cchi(1),cclo(2):cchi(2))
    integer         , intent(in   ) :: msk( mlo(1): mhi(1), mlo(2): mhi(2))

    integer, dimension(2) :: lo, hi, glo, ghi
    integer :: i, j, ii, jj, step
    real(amrex_real), parameter :: w1 = 9.d0/64.d0
    real(amrex_real), parameter :: w2 = 3.d0/64.d0
    real(amrex_real), parameter :: w3 = 1.d0/64.d0

    lo = 2*clo
    hi = 2*chi
    glo = 2*cglo
    ghi = 2*cghi

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
             rhs(i,j) = rhs(i,j) &
                  + w1*(cc(ii-1,jj-1)+cc(ii  ,jj-1)+cc(ii-1,jj  )+cc(ii  ,jj  )) &
                  + w2*(cc(ii-2,jj-1)+cc(ii+1,jj-1)+cc(ii-2,jj  )+cc(ii+1,jj  ) &
                  &    +cc(ii-1,jj-2)+cc(ii  ,jj-2)+cc(ii-1,jj+1)+cc(ii  ,jj+1)) &
                  + w3*(cc(ii-2,jj-2)+cc(ii+1,jj-2)+cc(ii-2,jj+1)+cc(ii+1,jj+1))
          end if
       end do
    end do

  end subroutine amrex_mlndlap_rhcc_fine_contrib


  subroutine amrex_mlndlap_rhcc_crse_contrib (lo, hi, crhs, rlo, rhi, rhcc, clo, chi, &
       dmsk, mlo, mhi, ndmsk, nmlo, nmhi, ccmsk, cmlo, cmhi) &
       bind(c,name='amrex_mlndlap_rhcc_crse_contrib')
    integer, dimension(2), intent(in) :: lo, hi, rlo, rhi, clo, chi, mlo, mhi, &
         nmlo, nmhi, cmlo, cmhi
    real(amrex_real), intent(inout) ::  crhs( rlo(1): rhi(1), rlo(2): rhi(2))
    real(amrex_real), intent(in   ) ::  rhcc( clo(1): chi(1), clo(2): chi(2))
    integer         , intent(in   ) ::  dmsk( mlo(1): mhi(1), mlo(2): mhi(2))
    integer         , intent(in   ) :: ndmsk(nmlo(1):nmhi(1),nmlo(2):nmhi(2))
    integer         , intent(in   ) :: ccmsk(cmlo(1):cmhi(1),cmlo(2):cmhi(2))

    integer :: i,j

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (dmsk(i,j) .ne. dirichlet) then
             if (ndmsk(i,j) .eq. crse_fine_node) then
                crhs(i,j) = crhs(i,j) + 0.25d0 * &
                     ( (1.d0-ccmsk(i-1,j-1)) * rhcc(i-1,j-1) &
                     + (1.d0-ccmsk(i  ,j-1)) * rhcc(i  ,j-1) &
                     + (1.d0-ccmsk(i-1,j  )) * rhcc(i-1,j  ) &
                     + (1.d0-ccmsk(i  ,j  )) * rhcc(i  ,j  ))
             end if
          end if
       end do
    end do

  end subroutine amrex_mlndlap_rhcc_crse_contrib


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


  subroutine amrex_mlndlap_stencil_rap (lo, hi, csten, clo, chi, fsten, flo, fhi) &
       bind(c,name='amrex_mlndlap_stencil_rap')
    integer, dimension(2), intent(in) :: lo, hi, clo, chi, flo, fhi
    real(amrex_real), intent(inout) :: csten(clo(1):chi(1),clo(2):chi(2),5)
    real(amrex_real), intent(in   ) :: fsten(flo(1):fhi(1),flo(2):fhi(2),5)
    
    integer :: i,j, ii, jj
    real(amrex_real) :: ap(-1:1,-1:1), p(-1:1,-1:1), cross1, cross2

    do    j = lo(2), hi(2)
       jj = 2*j
       do i = lo(1), hi(1)
          ii = 2*i

          ap = 0.d0
          p = 0.d0
          
          ! csten(i,j,2)
          p(-1,-1) = interp_from_pp_to(ii+1,jj-1)
          p( 0,-1) = interp_from_0p_to(ii+2,jj-1)
          p(-1, 0) = interp_from_p0_to(ii+1,jj  )
          p( 0, 0) = 1.d0
          p(-1, 1) = interp_from_pm_to(ii+1,jj+1)
          p( 0, 1) = interp_from_0m_to(ii+2,jj+1)
          
          ap(0,-1) = Ap0(ii,jj-1)*p(-1,-1) + App(ii,jj-1)*p(-1,0)
          ap(1,-1) = A00(ii+1,jj-1)*p(-1,-1) + Ap0(ii+1,jj-1)*p(0,-1) &
               + A0p(ii+1,jj-1)*p(-1,0) + App(ii+1,jj-1)*p(0,0)
          ap(0,0) = Apm(ii,jj)*p(-1,-1) + Ap0(ii,jj)*p(-1,0) + App(ii,jj)*p(-1,1)
          ap(1,0) = A0m(ii+1,jj)*p(-1,-1) + Apm(ii+1,jj)*p(0,-1) &
               + A00(ii+1,jj)*p(-1,0) + Ap0(ii+1,jj)*p(0,0) &
               + A0p(ii+1,jj)*p(-1,1) + App(ii+1,jj)*p(0,1)
          ap(0,1) = Apm(ii,jj+1)*p(-1,0) + Ap0(ii,jj+1)*p(-1,1)
          ap(1,1) = A0m(ii+1,jj+1)*p(-1,0) + Apm(ii+1,jj+1)*p(0,0) &
               + A00(ii+1,jj+1)*p(-1,1) + Ap0(ii+1,jj+1)*p(0,1)
          
          csten(i,j,2) = 0.25d0*(restrict_from_0m_to(ii,jj)*ap(0,-1) &
               + restrict_from_pm_to(ii,jj)*ap(1,-1) &
               + ap(0,0) &
               + restrict_from_p0_to(ii,jj)*ap(1,0) &
               + restrict_from_0p_to(ii,jj)*ap(0,1) &
               + restrict_from_pp_to(ii,jj)*ap(1,1))
          
          ! csten(i,j,3)
          p(-1,-1) = interp_from_pp_to(ii-1,jj+1)
          p( 0,-1) = interp_from_0p_to(ii  ,jj+1)
          p( 1,-1) = interp_from_mp_to(ii+1,jj+1)
          p(-1, 0) = interp_from_p0_to(ii-1,jj+2)
          p( 0, 0) = 1.d0
          p( 1, 0) = interp_from_m0_to(ii+1,jj+2)
          
          ap(-1,0) = A0p(ii-1,jj)*p(-1,-1) + App(ii-1,jj)*p(0,-1)
          ap(0,0) = Amp(ii,jj)*p(-1,-1) + A0p(ii,jj)*p(0,-1) + App(ii,jj)*p(1,-1)
          ap(1,0) = Amp(ii+1,jj)*p(0,-1) + A0p(ii+1,jj)*p(1,-1)
          ap(-1,1) = A00(ii-1,jj+1)*p(-1,-1) + Ap0(ii-1,jj+1)*p(0,-1) &
               + A0p(ii-1,jj+1)*p(-1,0) + App(ii-1,jj+1)*p(0,0)
          ap(0,1) = Am0(ii,jj+1)*p(-1,-1) + A00(ii,jj+1)*p(0,-1) + Ap0(ii,jj+1)*p(1,-1) &
               + Amp(ii,jj+1)*p(-1,0) + A0p(ii,jj+1)*p(0,0) + App(ii,jj+1)*p(1,0)
          ap(1,1) = Am0(ii+1,jj+1)*p(0,-1) + A00(ii+1,jj+1)*p(1,-1) &
               + Amp(ii+1,jj+1)*p(0,0) + A0p(ii+1,jj+1)*p(1,0)
          
          csten(i,j,3) = 0.25d0*(restrict_from_m0_to(ii,jj)*ap(-1,0) &
               + ap(0,0) &
               + restrict_from_p0_to(ii,jj)*ap(1,0) &
               + restrict_from_mp_to(ii,jj)*ap(-1,1) &
               + restrict_from_0p_to(ii,jj)*ap(0,1) &
               + restrict_from_pp_to(ii,jj)*ap(1,1))
          
          ! csten(i,j,4)
          p(-1,-1) = interp_from_pp_to(ii+1,jj+1)
          p( 0,-1) = interp_from_0p_to(ii+2,jj+1)
          p(-1, 0) = interp_from_p0_to(ii+1,jj+2)
          p( 0, 0) = 1.d0
          
          ap(0,0) = App(ii,jj)*p(-1,-1)
          ap(1,0) = A0p(ii+1,jj)*p(-1,-1) + App(ii+1,jj)*p(0,-1)
          ap(0,1) = Ap0(ii,jj+1)*p(-1,-1) + App(ii,jj+1)*p(-1,0)
          ap(1,1) = A00(ii+1,jj+1)*p(-1,-1) + Ap0(ii+1,jj+1)*p(0,-1) &
               + A0p(ii+1,jj+1)*p(-1,0) + App(ii+1,jj+1)*p(0,0)
          
          cross1 = 0.25d0*(ap(0,0) &
               + restrict_from_p0_to(ii,jj)*ap(1,0) &
               + restrict_from_0p_to(ii,jj)*ap(0,1) &
               + restrict_from_pp_to(ii,jj)*ap(1,1))

          p(0,-1) = interp_from_0p_to(ii,jj+1)
          p(1,-1) = interp_from_mp_to(ii+1,jj+1)
          p(0, 0) = 1.d0
          p(1, 0) = interp_from_m0_to(ii+1,jj+2)

          ap(-1,0) = Amp(ii+1,jj)*p(0,-1) + A0p(ii+1,jj)*p(1,-1)
          ap( 0,0) = Amp(ii+2,jj)*p(1,-1)
          ap(-1,1) = Am0(ii+1,jj+1)*p(0,-1) + A00(ii+1,jj+1)*p(1,-1) + Amp(ii+1,jj+1)*p(0,0) &
               + A0p(ii+1,jj+1)*p(1,0)
          ap( 0,1) = Am0(ii+2,jj+1)*p(1,-1) + Amp(ii+2,jj+1)*p(1,0)

          cross2 = 0.25*(ap(0,0) &
               + restrict_from_m0_to(ii+2,jj)*ap(-1,0) &
               + restrict_from_mp_to(ii+2,jj)*ap(-1,1) &
               + restrict_from_0p_to(ii+2,jj)*ap( 0,1))

          csten(i,j,4) = 0.5d0*(cross1+cross2)

       end do
    end do

  contains

    elemental function interp_from_mm_to (i,j) result(p)
      integer, intent(in) :: i,j
      real(amrex_real) :: p, wxm, wym, wmm
      wxm = abs(fsten(i-1,j  ,2))/(abs(fsten(i-1,j-1,4))+abs(fsten(i-1,j  ,4))+eps)
      wym = abs(fsten(i  ,j-1,3))/(abs(fsten(i-1,j-1,4))+abs(fsten(i  ,j-1,4))+eps)
      wmm = abs(fsten(i-1,j-1,4)) * (1.d0 + wxm + wym)
      p = wmm * fsten(i,j,5)
    end function interp_from_mm_to

    elemental function interp_from_mp_to (i,j) result(p)
      integer, intent(in) :: i,j
      real(amrex_real) :: p, wxm, wyp, wmp
      wxm = abs(fsten(i-1,j  ,2))/(abs(fsten(i-1,j-1,4))+abs(fsten(i-1,j  ,4))+eps)
      wyp = abs(fsten(i  ,j  ,3))/(abs(fsten(i-1,j  ,4))+abs(fsten(i  ,j  ,4))+eps)
      wmp = abs(fsten(i-1,j  ,4)) *(1.d0 + wxm + wyp)
      p = wmp * fsten(i,j,5)
    end function interp_from_mp_to

    elemental function interp_from_pm_to (i,j) result(p)
      integer, intent(in) :: i,j
      real(amrex_real) :: p, wxp, wym, wpm
      wxp = abs(fsten(i  ,j  ,2))/(abs(fsten(i  ,j-1,4))+abs(fsten(i  ,j  ,4))+eps)
      wym = abs(fsten(i  ,j-1,3))/(abs(fsten(i-1,j-1,4))+abs(fsten(i  ,j-1,4))+eps)
      wpm = abs(fsten(i  ,j-1,4)) * (1.d0 + wxp + wym)
      p = wpm * fsten(i,j,5)
    end function interp_from_pm_to

    elemental function interp_from_pp_to (i,j) result(p)
      integer, intent(in) :: i,j
      real(amrex_real) :: p, wxp, wyp, wpp
      wxp = abs(fsten(i  ,j  ,2))/(abs(fsten(i  ,j-1,4))+abs(fsten(i  ,j  ,4))+eps)
      wyp = abs(fsten(i  ,j  ,3))/(abs(fsten(i-1,j  ,4))+abs(fsten(i  ,j  ,4))+eps)
      wpp = abs(fsten(i  ,j  ,4)) * (1.d0 + wxp + wyp)
      p = wpp * fsten(i,j,5)
    end function interp_from_pp_to

    elemental function interp_from_m0_to (i,j) result(p)
      integer, intent(in) :: i,j
      real(amrex_real) :: p
      p = abs(fsten(i-1,j,2))/(abs(fsten(i-1,j,2))+abs(fsten(i,j,2))+eps)
    end function interp_from_m0_to

    elemental function interp_from_p0_to (i,j) result(p)
      integer, intent(in) :: i,j
      real(amrex_real) :: p
      p = abs(fsten(i,j,2))/(abs(fsten(i-1,j,2))+abs(fsten(i,j,2))+eps)
    end function interp_from_p0_to
    
    elemental function interp_from_0m_to (i,j) result(p)
      integer, intent(in) :: i,j
      real(amrex_real) :: p
      p = abs(fsten(i,j-1,3))/(abs(fsten(i,j-1,3))+abs(fsten(i,j,3))+eps)
    end function interp_from_0m_to

    elemental function interp_from_0p_to (i,j) result(p)
      integer, intent(in) :: i,j
      real(amrex_real) :: p
      p = abs(fsten(i,j,3))/(abs(fsten(i,j-1,3))+abs(fsten(i,j,3))+eps)
    end function interp_from_0p_to

    elemental real(amrex_real) function Amm (i,j)
      integer, intent(in) :: i,j
      Amm = fsten(i-1,j-1,4)
    end function Amm

    elemental real(amrex_real) function A0m (i,j)
      integer, intent(in) :: i,j
      A0m = fsten(i,j-1,3)
    end function A0m

    elemental real(amrex_real) function Apm (i,j)
      integer, intent(in) :: i,j
      Apm = fsten(i,j-1,4)
    end function Apm

    elemental real(amrex_real) function Am0 (i,j)
      integer, intent(in) :: i,j
      Am0 = fsten(i-1,j,2)
    end function Am0

    elemental real(amrex_real) function A00 (i,j)
      integer, intent(in) :: i,j
      A00 = fsten(i,j,1)
    end function A00

    elemental real(amrex_real) function Ap0 (i,j)
      integer, intent(in) :: i,j
      Ap0 = fsten(i,j,2)
    end function Ap0

    elemental real(amrex_real) function Amp (i,j)
      integer, intent(in) :: i,j
      Amp = fsten(i-1,j,4)
    end function Amp

    elemental real(amrex_real) function A0p (i,j)
      integer, intent(in) :: i,j
      A0p = fsten(i,j,3)
    end function A0p

    elemental real(amrex_real) function App (i,j)
      integer, intent(in) :: i,j
      App = fsten(i,j,4)
    end function App

    elemental function restrict_from_mm_to (ii,jj) result(r)
      integer, intent(in) :: ii,jj
      real(amrex_real) :: r, wxp, wyp, wpp
      wxp = abs(fsten(ii-1,jj-1,2))/(abs(fsten(ii-1,jj-2,4))+abs(fsten(ii-1,jj-1,4))+eps)
      wyp = abs(fsten(ii-1,jj-1,3))/(abs(fsten(ii-2,jj-1,4))+abs(fsten(ii-1,jj-1,4))+eps)
      wpp = abs(fsten(ii-1,jj-1,4))*(1.d0+wxp+wyp)
      r = wpp * fsten(ii-1,jj-1,5)
    end function restrict_from_mm_to

    elemental function restrict_from_0m_to (ii,jj) result(r)
      integer, intent(in) :: ii,jj
      real(amrex_real) :: r
      r = abs(fsten(ii,jj-1,3))/(abs(fsten(ii,jj-2,3))+abs(fsten(ii,jj-1,3))+eps)
    end function restrict_from_0m_to

    elemental function restrict_from_pm_to (ii,jj) result(r)
      integer, intent(in) :: ii,jj
      real(amrex_real) :: r, wxm, wyp, wmp
      wxm = abs(fsten(ii  ,jj-1,2))/(abs(fsten(ii,jj-2,4))+abs(fsten(ii  ,jj-1,4))+eps)
      wyp = abs(fsten(ii+1,jj-1,3))/(abs(fsten(ii,jj-1,4))+abs(fsten(ii+1,jj-1,4))+eps)
      wmp = abs(fsten(ii  ,jj-1,4)) *(1.d0 + wxm + wyp)
      r = wmp * fsten(ii+1,jj-1,5)
    end function restrict_from_pm_to

    elemental function restrict_from_m0_to (ii,jj) result(r)
      integer, intent(in) :: ii,jj
      real(amrex_real) :: r
      r = abs(fsten(ii-1,jj,2))/(abs(fsten(ii-2,jj,2))+abs(fsten(ii-1,jj,2))+eps)
    end function restrict_from_m0_to

    elemental function restrict_from_p0_to (ii,jj) result(r)
      integer, intent(in) :: ii,jj
      real(amrex_real) :: r
      r = abs(fsten(ii,jj,2))/(abs(fsten(ii,jj,2))+abs(fsten(ii+1,jj,2))+eps)
    end function restrict_from_p0_to

    elemental function restrict_from_mp_to (ii,jj) result(r)
      integer, intent(in) :: ii,jj
      real(amrex_real) :: r, wxp, wym, wpm
      wxp = abs(fsten(ii-1,jj+1,2))/(abs(fsten(ii-1,jj,4))+abs(fsten(ii-1,jj+1,4))+eps)
      wym = abs(fsten(ii-1,jj  ,3))/(abs(fsten(ii-2,jj,4))+abs(fsten(ii-1,jj  ,4))+eps)
      wpm = abs(fsten(ii-1,jj  ,4)) * (1.d0 + wxp + wym)
      r = wpm * fsten(ii-1,jj+1,5)
    end function restrict_from_mp_to

    elemental function restrict_from_0p_to (ii,jj) result(r)
      integer, intent(in) :: ii,jj
      real(amrex_real) :: r
      r = abs(fsten(ii,jj,3))/(abs(fsten(ii,jj,3))+abs(fsten(ii,jj+1,3))+eps)
    end function restrict_from_0p_to
    
    elemental function restrict_from_pp_to (ii,jj) result(r)
      integer, intent(in) :: ii,jj
      real(amrex_real) :: r, wxm, wym, wmm
      wxm = abs(fsten(ii  ,jj+1,2))/(abs(fsten(ii  ,jj  ,4))+abs(fsten(ii  ,jj+1,4))+eps)
      wym = abs(fsten(ii+1,jj  ,3))/(abs(fsten(ii  ,jj  ,4))+abs(fsten(ii+1,jj  ,4))+eps)
      wmm = abs(fsten(ii  ,jj  ,4)) * (1.d0 + wxm + wym)
      r = wmm * fsten(ii+1,jj+1,5)
    end function restrict_from_pp_to

  end subroutine amrex_mlndlap_stencil_rap


#ifdef AMREX_USE_EB

  subroutine amrex_mlndlap_set_integral (lo, hi, intg, glo, ghi) &
       bind(c,name='amrex_mlndlap_set_integral')
    integer, dimension(2) :: lo, hi, glo, ghi
    real(amrex_real), intent(inout) :: intg(glo(1):ghi(1),glo(2):ghi(2),n_Sintg)
    integer :: i,j
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          intg(i,j,i_S_x ) = zero
          intg(i,j,i_S_y ) = zero
          intg(i,j,i_S_x2) = twelfth
          intg(i,j,i_S_y2) = twelfth
          intg(i,j,i_S_xy) = zero
       end do
    end do
  end subroutine amrex_mlndlap_set_integral

  subroutine amrex_mlndlap_set_integral_eb (lo, hi, intg, glo, ghi, flag, flo, fhi, &
       vol, vlo, vhi, ax, axlo, axhi, ay, aylo, ayhi, bcen, blo, bhi) &
       bind(c,name='amrex_mlndlap_set_integral_eb')
    use amrex_ebcellflag_module, only : is_single_valued_cell, is_regular_cell, is_covered_cell
    integer, dimension(2) :: lo, hi, glo, ghi, flo, fhi, axlo, vlo, vhi, axhi, aylo, ayhi, blo, bhi
    real(amrex_real), intent(inout) :: intg( glo(1): ghi(1), glo(2): ghi(2),n_Sintg)
    real(amrex_real), intent(in   ) :: vol ( vlo(1): vhi(1), vlo(2): vhi(2))
    real(amrex_real), intent(in   ) :: ax  (axlo(1):axhi(1),axlo(2):axhi(2))
    real(amrex_real), intent(in   ) :: ay  (aylo(1):ayhi(1),aylo(2):ayhi(2))
    real(amrex_real), intent(in   ) :: bcen( blo(1): bhi(1), blo(2): bhi(2),2)
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2))

    integer :: i,j
    real(amrex_real) :: Sx, Sx2, Sy, Sy2, Sxy ! integral of x, x2, y, y2 and xy
    real(amrex_real) :: axm, axp, aym, ayp, apnorm, apnorminv, anrmx, anrmy, bcx, bcy, kk, bb
    real(amrex_real) :: xmin, xmax, ymin, ymax
    real(amrex_real), parameter :: almostone = 1.d0 - 1.d2*epsilon(1._amrex_real)
    real(amrex_real), parameter :: sixteenth = 1.d0/16.d0, twentyfourth = 1.d0/24.d0

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (is_covered_cell(flag(i,j))) then

             intg(i,j,:) = zero

          else if (is_regular_cell(flag(i,j)) .or. vol(i,j).ge.almostone) then

             intg(i,j,i_S_x ) = zero
             intg(i,j,i_S_y ) = zero
             intg(i,j,i_S_x2) = twelfth
             intg(i,j,i_S_y2) = twelfth
             intg(i,j,i_S_xy) = zero

          else

             axm = ax(i,j)
             axp = ax(i+1,j)
             aym = ay(i,j)
             ayp = ay(i,j+1)

             apnorm = sqrt((axm-axp)**2 + (aym-ayp)**2)
             if (apnorm .eq. zero) then
                call amrex_abort("amrex_mlndlap_set_integral: we are in trouble")
             end if

             apnorminv = 1.d0/apnorm
             anrmx = (axm-axp) * apnorminv  ! pointing to the wall
             anrmy = (aym-ayp) * apnorminv

             bcx = bcen(i,j,1)
             bcy = bcen(i,j,2)

             if (anrmx .eq. zero) then
                Sx = zero
                Sx2 = twentyfourth*(axm+axp)
                Sxy = zero
             else if (anrmy .eq. zero) then
                Sx  = eighth     *(axp-axm) + anrmx*half*bcx**2
                Sx2 = twentyfourth*(axp+axm) + anrmx*third*bcx**3
                Sxy = zero
             else
                if (anrmx .gt. zero) then
                   xmin = -half + min(aym,ayp)
                   xmax = -half + max(aym,ayp)
                else
                   xmin = half - max(aym,ayp)
                   xmax = half - min(aym,ayp)
                end if
                Sx  = eighth     *(axp-axm) + (anrmx/abs(anrmy))*sixth  *(xmax**3-xmin**3)
                Sx2 = twentyfourth*(axp+axm) + (anrmx/abs(anrmy))*twelfth*(xmax**4-xmin**4)

                kk = -anrmx/anrmy
                bb = bcy-kk*bcx
                Sxy = eighth*kk*kk*(xmax**4-xmin**4) + third*kk*bb*(xmax**3-xmin**3) &
                     + (fourth*bb*bb-sixteenth)*(xmax**2-xmin**2)
                sxy = sxy * sign(one,anrmy)
             end if

             if (anrmy .eq. zero) then
                Sy = zero
                Sy2 = twentyfourth*(aym+ayp)
             else if (anrmx .eq. zero) then
                Sy  = eighth     *(ayp-aym) + anrmy*half*bcy**2
                Sy2 = twentyfourth*(ayp+aym) + anrmy*third*bcy**3
             else
                if (anrmy .gt. zero) then
                   ymin = -half + min(axm,axp)
                   ymax = -half + max(axm,axp)
                else
                   ymin = half - max(axm,axp)
                   ymax = half - min(axm,axp)
                end if
                Sy  = eighth     *(ayp-aym) + (anrmy/abs(anrmx))*sixth  *(ymax**3-ymin**3)
                Sy2 = twentyfourth*(ayp+aym) + (anrmy/abs(anrmx))*twelfth*(ymax**4-ymin**4)
             end if

             intg(i,j,i_S_x ) = Sx
             intg(i,j,i_S_y ) = Sy
             intg(i,j,i_S_x2) = Sx2
             intg(i,j,i_S_y2) = Sy2
             intg(i,j,i_S_xy) = Sxy

          end if
       end do
    end do
  end subroutine amrex_mlndlap_set_integral_eb


  subroutine amrex_mlndlap_set_connection (lo, hi, conn, clo, chi, intg, glo, ghi, flag, flo, fhi, &
       vol, vlo, vhi) bind(c,name='amrex_mlndlap_set_connection')
    use amrex_ebcellflag_module, only : is_single_valued_cell, is_regular_cell, is_covered_cell
    integer, dimension(2), intent(in) :: lo, hi, clo, chi, glo, ghi, flo, fhi, vlo, vhi
    real(amrex_real), intent(inout) :: conn( clo(1): chi(1), clo(2): chi(2),6)
    real(amrex_real), intent(in   ) :: intg( glo(1): ghi(1), glo(2): ghi(2),n_Sintg)
    real(amrex_real), intent(in   ) :: vol ( vlo(1): vhi(1), vlo(2): vhi(2))
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2))

    integer :: i,j
    real(amrex_real), parameter :: almostone = 1.d0 - 1.d2*epsilon(1._amrex_real)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (is_covered_cell(flag(i,j))) then

             conn(i,j,:) = zero

          else if (is_regular_cell(flag(i,j)) .or. vol(i,j).ge.almostone) then

             conn(i,j,:) = 1.d0

          else

             ! Note that these are normalized so that they equal 1 in the case of a regular cell

             conn(i,j,1) = 3.d0*(.25d0*vol(i,j) + intg(i,j,i_S_y2) - intg(i,j,i_S_y))
             conn(i,j,2) = 6.d0*(.25d0*vol(i,j) - intg(i,j,i_S_y2))
             conn(i,j,3) = 3.d0*(.25d0*vol(i,j) + intg(i,j,i_S_y2) + intg(i,j,i_S_y))

             conn(i,j,4) = 3.d0*(.25d0*vol(i,j) + intg(i,j,i_S_x2) - intg(i,j,i_S_x))
             conn(i,j,5) = 6.d0*(.25d0*vol(i,j) - intg(i,j,i_S_x2))
             conn(i,j,6) = 3.d0*(.25d0*vol(i,j) + intg(i,j,i_S_x2) + intg(i,j,i_S_x))

          end if
       end do
    end do
  end subroutine amrex_mlndlap_set_connection


  subroutine amrex_mlndlap_set_stencil_eb (lo, hi, sten, tlo, thi, sigma, glo, ghi, &
       conn, clo, chi, dxinv) bind(c,name='amrex_mlndlap_set_stencil_eb')
    integer, dimension(2), intent(in) :: lo, hi, tlo, thi, glo, ghi, clo, chi
    real(amrex_real), intent(inout) ::  sten(tlo(1):thi(1),tlo(2):thi(2),5)
    real(amrex_real), intent(in   ) :: sigma(glo(1):ghi(1),glo(2):ghi(2))
    real(amrex_real), intent(in   ) ::  conn(clo(1):chi(1),clo(2):chi(2),6)
    real(amrex_real), intent(in) :: dxinv(2)

    integer :: i, j
    real(amrex_real) :: facx, facy

    facx = (1.d0/6.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/6.d0)*dxinv(2)*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          sten(i,j,2) = 2.d0*facx*(sigma(i,j-1)*conn(i,j-1,3)+sigma(i,j)*conn(i,j,1)) &
               &            -facy*(sigma(i,j-1)*conn(i,j-1,5)+sigma(i,j)*conn(i,j,5))
          sten(i,j,3) = 2.d0*facy*(sigma(i-1,j)*conn(i-1,j,6)+sigma(i,j)*conn(i,j,4)) &
               &            -facx*(sigma(i-1,j)*conn(i-1,j,2)+sigma(i,j)*conn(i,j,2))
          sten(i,j,4) = (facx*conn(i,j,2)+facy*conn(i,j,5))*sigma(i,j)
       end do
    end do

  end subroutine amrex_mlndlap_set_stencil_eb


  subroutine amrex_mlndlap_divu_eb (lo, hi, rhs, rlo, rhi, vel, vlo, vhi, vfrac, flo, fhi, &
       intg, glo, ghi, msk, mlo, mhi, dxinv) &
       bind(c,name='amrex_mlndlap_divu_eb')
    integer, dimension(2), intent(in) :: lo, hi, rlo, rhi, vlo, vhi, flo, fhi, glo, ghi, mlo, mhi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1),vlo(2):vhi(2),2)
    real(amrex_real), intent(in   ) :: vfrac(flo(1):fhi(1),flo(2):fhi(2))
    real(amrex_real), intent(in   ) :: intg(glo(1):ghi(1),glo(2):ghi(2),n_Sintg)
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))

    integer :: i,j
    real(amrex_real) :: facx, facy

    facx = half*dxinv(1)
    facy = half*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (msk(i,j) .ne. dirichlet) then
             rhs(i,j) = facx*(-vel(i-1,j-1,1)*(vfrac(i-1,j-1)+2.d0*intg(i-1,j-1,2)) &
                  &           +vel(i  ,j-1,1)*(vfrac(i  ,j-1)+2.d0*intg(i  ,j-1,2)) &
                  &           -vel(i-1,j  ,1)*(vfrac(i-1,j  )-2.d0*intg(i-1,j  ,2)) &
                  &           +vel(i  ,j  ,1)*(vfrac(i  ,j  )-2.d0*intg(i  ,j  ,2))) &
                  &   + facy*(-vel(i-1,j-1,2)*(vfrac(i-1,j-1)+2.d0*intg(i-1,j-1,1)) &
                  &           -vel(i  ,j-1,2)*(vfrac(i  ,j-1)-2.d0*intg(i  ,j-1,1)) &
                  &           +vel(i-1,j  ,2)*(vfrac(i-1,j  )+2.d0*intg(i-1,j  ,1)) &
                  &           +vel(i  ,j  ,2)*(vfrac(i  ,j  )-2.d0*intg(i  ,j  ,1)))
          else
             rhs(i,j) = zero
          end if
       end do
    end do

  end subroutine amrex_mlndlap_divu_eb


  subroutine amrex_mlndlap_mknewu_eb (lo, hi, u, ulo, uhi, p, plo, phi, sig, slo, shi, &
       vfrac, vlo, vhi, intg, glo, ghi, dxinv) bind(c,name='amrex_mlndlap_mknewu_eb')
    integer, dimension(2), intent(in) :: lo, hi, ulo, uhi, plo, phi, slo, shi, vlo, vhi, glo, ghi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) ::   u(ulo(1):uhi(1),ulo(2):uhi(2),2)
    real(amrex_real), intent(in   ) ::   p(plo(1):phi(1),plo(2):phi(2))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1),slo(2):shi(2))
    real(amrex_real), intent(in   )::vfrac(vlo(1):vhi(1),vlo(2):vhi(2))
    real(amrex_real), intent(in   ) ::intg(glo(1):ghi(1),glo(2):ghi(2),n_Sintg)

    integer :: i, j
    real(amrex_real) :: dpdx, dpdy, facx, facy, dpp

    facx = half*dxinv(1)
    facy = half*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (vfrac(i,j) .eq. zero) then
             u(i,j,1) = zero
             u(i,j,2) = zero
          else
             dpdx = facx*(-p(i,j)+p(i+1,j)-p(i,j+1)+p(i+1,j+1))
             dpdy = facy*(-p(i,j)-p(i+1,j)+p(i,j+1)+p(i+1,j+1))
             dpp = (p(i,j)+p(i+1,j+1)-p(i+1,j)-p(i,j+1))/vfrac(i,j)
             u(i,j,1) = u(i,j,1) - sig(i,j)*(dpdx + dxinv(1)*intg(i,j,2)*dpp)
             u(i,j,2) = u(i,j,2) - sig(i,j)*(dpdy + dxinv(2)*intg(i,j,1)*dpp)
          end if
       end do
    end do
  end subroutine amrex_mlndlap_mknewu_eb

  subroutine amrex_mlndlap_rhcc_eb (lo, hi, rhs, rlo, rhi, rhcc, clo, chi, vfrac, flo,fhi,&
       intg, glo, ghi, msk, mlo, mhi) &
       bind(c,name='amrex_mlndlap_rhcc_eb')
    integer, dimension(2) :: lo, hi, rlo, rhi, clo, chi, flo, fhi, glo, ghi, mlo, mhi
    real(amrex_real), intent(inout) :: rhs (rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(in   ) :: rhcc(clo(1):chi(1),clo(2):chi(2))
    real(amrex_real), intent(in   ) :: vfrac(flo(1):fhi(1),flo(2):fhi(2))
    real(amrex_real), intent(in   ) :: intg(glo(1):ghi(1),glo(2):ghi(2),n_Sintg)
    integer,          intent(in   ) :: msk (mlo(1):mhi(1),mlo(2):mhi(2))

    integer :: i,j

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (msk(i,j) .ne. dirichlet) then
             rhs(i,j) = &
                  rhcc(i  ,j  )*(0.25d0*vfrac(i  ,j  )-intg(i  ,j  ,i_S_x)-intg(i  ,j  ,i_S_y)+intg(i  ,j  ,i_S_xy)) + &
                  rhcc(i-1,j  )*(0.25d0*vfrac(i-1,j  )+intg(i-1,j  ,i_S_x)-intg(i-1,j  ,i_S_y)-intg(i-1,j  ,i_S_xy)) + &
                  rhcc(i-1,j-1)*(0.25d0*vfrac(i-1,j-1)+intg(i-1,j-1,i_S_x)+intg(i-1,j-1,i_S_y)+intg(i-1,j-1,i_S_xy)) + &
                  rhcc(i  ,j-1)*(0.25d0*vfrac(i  ,j-1)-intg(i  ,j-1,i_S_x)+intg(i  ,j-1,i_S_y)-intg(i  ,j-1,i_S_xy))
          else
             rhs(i,j) = 0.d0
          end if
       end do
    end do
  end subroutine amrex_mlndlap_rhcc_eb

#endif

end module amrex_mlnodelap_2d_module
