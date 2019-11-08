
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

  real(amrex_real), private, parameter :: eps = 1.d-100

  integer, private, parameter :: ist_000 = 1
  integer, private, parameter :: ist_p00 = 2
  integer, private, parameter :: ist_0p0 = 3
  integer, private, parameter :: ist_00p = 4
  integer, private, parameter :: ist_pp0 = 5
  integer, private, parameter :: ist_p0p = 6
  integer, private, parameter :: ist_0pp = 7
  integer, private, parameter :: ist_ppp = 8
  integer, private, parameter :: ist_inv = 9
  integer, private, parameter :: n_sten = 9

#ifdef AMREX_USE_EB
  integer, private, parameter :: i_S_x     = 1
  integer, private, parameter :: i_S_y     = 2
  integer, private, parameter :: i_S_z     = 3
  integer, private, parameter :: i_S_x2    = 4
  integer, private, parameter :: i_S_y2    = 5
  integer, private, parameter :: i_S_z2    = 6
  integer, private, parameter :: i_S_x_y   = 7
  integer, private, parameter :: i_S_x_z   = 8
  integer, private, parameter :: i_S_y_z   = 9
  integer, private, parameter :: i_S_x2_y  = 10
  integer, private, parameter :: i_S_x2_z  = 11
  integer, private, parameter :: i_S_x_y2  = 12
  integer, private, parameter :: i_S_y2_z  = 13
  integer, private, parameter :: i_S_x_z2  = 14
  integer, private, parameter :: i_S_y_z2  = 15
  integer, private, parameter :: i_S_x2_y2 = 16
  integer, private, parameter :: i_S_x2_z2 = 17
  integer, private, parameter :: i_S_y2_z2 = 18
  integer, private, parameter :: i_S_xyz   = 19
  integer, private, parameter :: n_Sintg   = 19

  integer, private, parameter :: i_c_xmym = 1
  integer, private, parameter :: i_c_xmyb = 2
  integer, private, parameter :: i_c_xmyp = 3
  integer, private, parameter :: i_c_xbym = 4
  integer, private, parameter :: i_c_xbyb = 5
  integer, private, parameter :: i_c_xbyp = 6
  integer, private, parameter :: i_c_xpym = 7
  integer, private, parameter :: i_c_xpyb = 8
  integer, private, parameter :: i_c_xpyp = 9
  integer, private, parameter :: i_c_xmzm = 10
  integer, private, parameter :: i_c_xmzb = 11
  integer, private, parameter :: i_c_xmzp = 12
  integer, private, parameter :: i_c_xbzm = 13
  integer, private, parameter :: i_c_xbzb = 14
  integer, private, parameter :: i_c_xbzp = 15
  integer, private, parameter :: i_c_xpzm = 16
  integer, private, parameter :: i_c_xpzb = 17
  integer, private, parameter :: i_c_xpzp = 18
  integer, private, parameter :: i_c_ymzm = 19
  integer, private, parameter :: i_c_ymzb = 20
  integer, private, parameter :: i_c_ymzp = 21
  integer, private, parameter :: i_c_ybzm = 22
  integer, private, parameter :: i_c_ybzb = 23
  integer, private, parameter :: i_c_ybzp = 24
  integer, private, parameter :: i_c_ypzm = 25
  integer, private, parameter :: i_c_ypzb = 26
  integer, private, parameter :: i_c_ypzp = 27
  integer, private, parameter :: n_conn = 27
#endif

  private
  public :: &
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

#ifdef AMREX_USE_EB
  public:: amrex_mlndlap_mknewu_eb, amrex_mlndlap_rhcc_eb
#endif

contains

  function amrex_mlndlap_any_fine_sync_cells (lo, hi, msk, mlo, mhi, fine_flag) result(r) &
       bind(c,name='amrex_mlndlap_any_fine_sync_cells')
    integer :: r
    integer, dimension(3), intent(in) :: lo, hi, mlo, mhi
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    integer, intent(in) :: fine_flag

    integer :: i,j,k

    r = 0

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .eq. fine_flag) then
                r = 1
                goto 100
             end if
          end do
       end do
    end do
100 continue
  end function amrex_mlndlap_any_fine_sync_cells


  subroutine amrex_mlndlap_divu_fine_contrib (clo, chi, cglo, cghi, rhs, rlo, rhi, &
       vel, vlo, vhi, frh, flo, fhi, msk, mlo, mhi, dxinv) &
       bind(c,name='amrex_mlndlap_divu_fine_contrib')
    integer, dimension(3), intent(in) :: clo, chi, cglo, cghi, rlo, rhi, vlo, vhi, &
         flo, fhi, mlo, mhi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3)
    real(amrex_real), intent(inout) :: frh(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    integer         , intent(in   ) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer, dimension(3) :: lo, hi, glo, ghi, gtlo, gthi
    integer :: i, j, k, ii, jj, kk, step
    real(amrex_real) :: facx, facy, facz
    real(amrex_real), parameter :: rfd = 0.125d0
    real(amrex_real), parameter :: chip = 0.5d0
    real(amrex_real), parameter :: chip2 = 0.25d0
    real(amrex_real), parameter :: chip3 = 0.125d0

    ! note that dxinv is fine dxinv
    facx = 0.25d0*dxinv(1)
    facy = 0.25d0*dxinv(2)
    facz = 0.25d0*dxinv(3)

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
                frh(ii,jj,kk) = facx*(-vel(ii-1,jj-1,kk-1,1)+vel(ii,jj-1,kk-1,1) &
                     &                -vel(ii-1,jj  ,kk-1,1)+vel(ii,jj  ,kk-1,1) &
                     &                -vel(ii-1,jj-1,kk  ,1)+vel(ii,jj-1,kk  ,1) &
                     &                -vel(ii-1,jj  ,kk  ,1)+vel(ii,jj  ,kk  ,1)) &
                     &        + facy*(-vel(ii-1,jj-1,kk-1,2)-vel(ii,jj-1,kk-1,2) &
                     &                +vel(ii-1,jj  ,kk-1,2)+vel(ii,jj  ,kk-1,2) &
                     &                -vel(ii-1,jj-1,kk  ,2)-vel(ii,jj-1,kk  ,2) &
                     &                +vel(ii-1,jj  ,kk  ,2)+vel(ii,jj  ,kk  ,2)) &
                     &        + facz*(-vel(ii-1,jj-1,kk-1,3)-vel(ii,jj-1,kk-1,3) &
                     &                -vel(ii-1,jj  ,kk-1,3)-vel(ii,jj  ,kk-1,3) &
                     &                +vel(ii-1,jj-1,kk  ,3)+vel(ii,jj-1,kk  ,3) &
                     &                +vel(ii-1,jj  ,kk  ,3)+vel(ii,jj  ,kk  ,3))
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
                rhs(i,j,k) = rhs(i,j,k) + rfd*(frh(ii,jj,kk) &
                     + chip*(frh(ii,jj,kk-1)+frh(ii,jj,kk+1) &
                     &      +frh(ii,jj-1,kk)+frh(ii,jj+1,kk) &
                     &      +frh(ii-1,jj,kk)+frh(ii+1,jj,kk)) &
                     + chip2*(frh(ii,jj-1,kk-1)+frh(ii,jj+1,kk-1)+frh(ii,jj-1,kk+1)+frh(ii,jj+1,kk+1) &
                     &       +frh(ii-1,jj,kk-1)+frh(ii+1,jj,kk-1)+frh(ii-1,jj,kk+1)+frh(ii+1,jj,kk+1) &
                     &       +frh(ii-1,jj-1,kk)+frh(ii+1,jj-1,kk)+frh(ii-1,jj+1,kk)+frh(ii+1,jj+1,kk)) &
                     + chip3*(frh(ii-1,jj-1,kk-1)+frh(ii+1,jj-1,kk-1) &
                     &       +frh(ii-1,jj+1,kk-1)+frh(ii+1,jj+1,kk-1) &
                     &       +frh(ii-1,jj-1,kk+1)+frh(ii+1,jj-1,kk+1) &
                     &       +frh(ii-1,jj+1,kk+1)+frh(ii+1,jj+1,kk+1)))
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_divu_fine_contrib

  subroutine amrex_mlndlap_divu_cf_contrib (lo, hi,  rhs, rlo, rhi, vel, vlo, vhi, dmsk, mlo, mhi, &
       ndmsk, nmlo, nmhi, ccmsk, cmlo, cmhi, fc, clo, chi, dxinv, ndlo, ndhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_divu_cf_contrib')
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, vlo, vhi, mlo, mhi, &
         nmlo, nmhi, cmlo, cmhi, clo, chi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: rhs  ( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) :: vel  ( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3),3)
    real(amrex_real), intent(in   ) :: fc   ( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3))
    integer         , intent(in   ) :: dmsk ( mlo(1): mhi(1), mlo(2): mhi(2), mlo(3): mhi(3))
    integer         , intent(in   ) :: ndmsk(nmlo(1):nmhi(1),nmlo(2):nmhi(2),nmlo(3):nmhi(3))
    integer         , intent(in   ) :: ccmsk(cmlo(1):cmhi(1),cmlo(2):cmhi(2),cmlo(3):cmhi(3))

    integer :: i,j,k
    real(amrex_real) :: facx, facy, facz

    facx = 0.25d0*dxinv(1)
    facy = 0.25d0*dxinv(2)
    facz = 0.25d0*dxinv(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (dmsk(i,j,k) .ne. dirichlet) then
                if (ndmsk(i,j,k) .eq. crse_fine_node) then
                   rhs(i,j,k) = fc(i,j,k)
                   if (ccmsk(i-1,j-1,k-1) .eq. crse_cell) then
                      rhs(i,j,k) = rhs(i,j,k) - facx*vel(i-1,j-1,k-1,1) &
                           &                  - facy*vel(i-1,j-1,k-1,2) &
                           &                  - facz*vel(i-1,j-1,k-1,3)
                   end if
                   if (ccmsk(i,j-1,k-1) .eq. crse_cell) then
                      rhs(i,j,k) = rhs(i,j,k) + facx*vel(i  ,j-1,k-1,1) &
                           &                  - facy*vel(i  ,j-1,k-1,2) &
                           &                  - facz*vel(i  ,j-1,k-1,3)
                   end if
                   if (ccmsk(i-1,j,k-1) .eq. crse_cell) then
                      rhs(i,j,k) = rhs(i,j,k) - facx*vel(i-1,j  ,k-1,1) &
                           &                  + facy*vel(i-1,j  ,k-1,2) &
                           &                  - facz*vel(i-1,j  ,k-1,3)
                   end if
                   if (ccmsk(i,j,k-1) .eq. crse_cell) then
                      rhs(i,j,k) = rhs(i,j,k) + facx*vel(i  ,j  ,k-1,1) &
                           &                  + facy*vel(i  ,j  ,k-1,2) &
                           &                  - facz*vel(i  ,j  ,k-1,3)
                   end if
                   if (ccmsk(i-1,j-1,k) .eq. crse_cell) then
                      rhs(i,j,k) = rhs(i,j,k) - facx*vel(i-1,j-1,k  ,1) &
                           &                  - facy*vel(i-1,j-1,k  ,2) &
                        &                     + facz*vel(i-1,j-1,k  ,3)
                   end if
                   if (ccmsk(i,j-1,k) .eq. crse_cell) then
                      rhs(i,j,k) = rhs(i,j,k) + facx*vel(i  ,j-1,k  ,1) &
                           &                  - facy*vel(i  ,j-1,k  ,2) &
                           &                  + facz*vel(i  ,j-1,k  ,3)
                   end if
                   if (ccmsk(i-1,j,k) .eq. crse_cell) then
                      rhs(i,j,k) = rhs(i,j,k) - facx*vel(i-1,j  ,k  ,1) &
                           &                  + facy*vel(i-1,j  ,k  ,2) &
                           &                  + facz*vel(i-1,j  ,k  ,3)
                   end if
                   if (ccmsk(i,j,k) .eq. crse_cell) then
                      rhs(i,j,k) = rhs(i,j,k) + facx*vel(i  ,j  ,k  ,1) &
                           &                  + facy*vel(i  ,j  ,k  ,2) &
                           &                  + facz*vel(i  ,j  ,k  ,3)
                   end if

                   if (i .eq. ndlo(1) .and. &
                        (    bclo(1) .eq. amrex_lo_neumann &
                        .or. bclo(1) .eq. amrex_lo_inflow)) then
                      rhs(i,j,k) = 2.d0*rhs(i,j,k)
                   else if (i.eq. ndhi(1) .and. &
                        (    bchi(1) .eq. amrex_lo_neumann &
                        .or. bchi(1) .eq. amrex_lo_inflow)) then
                      rhs(i,j,k) = 2.d0*rhs(i,j,k)
                   end if
                   
                   if (j .eq. ndlo(2) .and. &
                        (    bclo(2) .eq. amrex_lo_neumann &
                        .or. bclo(2) .eq. amrex_lo_inflow)) then
                      rhs(i,j,k) = 2.d0*rhs(i,j,k)                   
                   else if (j .eq. ndhi(2) .and. &
                        (    bchi(2) .eq. amrex_lo_neumann &
                        .or. bchi(2) .eq. amrex_lo_inflow)) then
                      rhs(i,j,k) = 2.d0*rhs(i,j,k)
                   end if

                   if (k .eq. ndlo(3) .and. &
                        (    bclo(3) .eq. amrex_lo_neumann &
                        .or. bclo(3) .eq. amrex_lo_inflow)) then
                      rhs(i,j,k) = 2.d0*rhs(i,j,k)                   
                   else if (k .eq. ndhi(3) .and. &
                        (    bchi(3) .eq. amrex_lo_neumann &
                        .or. bchi(3) .eq. amrex_lo_inflow)) then
                      rhs(i,j,k) = 2.d0*rhs(i,j,k)
                   end if

                end if
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_divu_cf_contrib


  subroutine amrex_mlndlap_rhcc_fine_contrib (clo, chi, cglo, cghi, rhs, rlo, rhi, &
       cc, cclo, cchi, msk, mlo, mhi) bind(c,name='amrex_mlndlap_rhcc_fine_contrib')
    integer, dimension(3), intent(in) :: clo, chi, cglo, cghi, rlo, rhi, cclo, cchi, mlo, mhi
    real(amrex_real), intent(inout) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) :: cc (cclo(1):cchi(1),cclo(2):cchi(2),cclo(3):cchi(3))
    integer         , intent(in   ) :: msk( mlo(1): mhi(1), mlo(2): mhi(2), mlo(3): mhi(3))

    integer, dimension(3) :: lo, hi, glo, ghi
    integer :: i, j, k, ii, jj, kk, step, ioff, joff, koff
    real(amrex_real), parameter :: fac(-2:1) = [1.d0/8.d0, 3.d0/8.d0, 3.d0/8.d0, 1.d0/8.d0]

    lo = 2*clo
    hi = 2*chi
    glo = 2*cglo
    ghi = 2*cghi

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
                do koff = -2, 1
                   do joff = -2, 1
                      do ioff = -2, 1
                         rhs(i,j,k) = rhs(i,j,k) + cc(ii+ioff,jj+joff,kk+koff) &
                              * fac(ioff)*fac(joff)*fac(koff)
                      end do
                   end do
                end do
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_rhcc_fine_contrib


  subroutine amrex_mlndlap_rhcc_crse_contrib (lo, hi, crhs, rlo, rhi, rhcc, clo, chi, &
       dmsk, mlo, mhi, ndmsk, nmlo, nmhi, ccmsk, cmlo, cmhi) &
       bind(c,name='amrex_mlndlap_rhcc_crse_contrib')
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, clo, chi, mlo, mhi, &
         nmlo, nmhi, cmlo, cmhi
    real(amrex_real), intent(inout) ::  crhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) ::  rhcc( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3))
    integer         , intent(in   ) ::  dmsk( mlo(1): mhi(1), mlo(2): mhi(2), mlo(3): mhi(3))
    integer         , intent(in   ) :: ndmsk(nmlo(1):nmhi(1),nmlo(2):nmhi(2),nmlo(3):nmhi(3))
    integer         , intent(in   ) :: ccmsk(cmlo(1):cmhi(1),cmlo(2):cmhi(2),cmlo(3):cmhi(3))

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (dmsk(i,j,k) .ne. dirichlet) then
                if (ndmsk(i,j,k) .eq. crse_fine_node) then
                   crhs(i,j,k) = crhs(i,j,k) + 0.125d0 * &
                        ( (1.d0-ccmsk(i-1,j-1,k-1)) * rhcc(i-1,j-1,k-1) &
                        + (1.d0-ccmsk(i  ,j-1,k-1)) * rhcc(i  ,j-1,k-1) &
                        + (1.d0-ccmsk(i-1,j  ,k-1)) * rhcc(i-1,j  ,k-1) &
                        + (1.d0-ccmsk(i  ,j  ,k-1)) * rhcc(i  ,j  ,k-1) &
                        + (1.d0-ccmsk(i-1,j-1,k  )) * rhcc(i-1,j-1,k  ) &
                        + (1.d0-ccmsk(i  ,j-1,k  )) * rhcc(i  ,j-1,k  ) &
                        + (1.d0-ccmsk(i-1,j  ,k  )) * rhcc(i-1,j  ,k  ) &
                        + (1.d0-ccmsk(i  ,j  ,k  )) * rhcc(i  ,j  ,k  ) )
                end if
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_rhcc_crse_contrib


  subroutine amrex_mlndlap_crse_resid (lo, hi, resid, rslo, rshi, rhs, rhlo, rhhi, msk, mlo, mhi, &
       ndlo, ndhi, bclo, bchi) bind(c, name='amrex_mlndlap_crse_resid')
    integer, dimension(3), intent(in) :: lo, hi, rslo, rshi, rhlo, rhhi, mlo, mhi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(inout) :: resid(rslo(1):rshi(1),rslo(2):rshi(2),rslo(3):rshi(3))
    real(amrex_real), intent(in   ) :: rhs  (rhlo(1):rhhi(1),rhlo(2):rhhi(2),rhlo(3):rhhi(3))
    integer         , intent(in   ) :: msk  ( mlo(1): mhi(1), mlo(2): mhi(2), mlo(3): mhi(3))

    integer :: i,j,k
    real(amrex_real) :: fac

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (any(msk(i-1:i,j-1:j,k-1:k).eq.0) .and. any(msk(i-1:i,j-1:j,k-1:k).eq.1)) then

                fac = 1.d0

                if (i .eq. ndlo(1) .and. &
                     (    bclo(1) .eq. amrex_lo_neumann &
                     .or. bclo(1) .eq. amrex_lo_inflow)) then
                   fac = 2.d0*fac
                else if (i.eq. ndhi(1) .and. &
                     (    bchi(1) .eq. amrex_lo_neumann &
                     .or. bchi(1) .eq. amrex_lo_inflow)) then
                   fac = 2.d0*fac
                end if
                
                if (j .eq. ndlo(2) .and. &
                     (    bclo(2) .eq. amrex_lo_neumann &
                     .or. bclo(2) .eq. amrex_lo_inflow)) then
                   fac = 2.d0*fac                   
                else if (j .eq. ndhi(2) .and. &
                     (    bchi(2) .eq. amrex_lo_neumann &
                     .or. bchi(2) .eq. amrex_lo_inflow)) then
                   fac = 2.d0*fac
                end if
                
                if (k .eq. ndlo(3) .and. &
                     (    bclo(3) .eq. amrex_lo_neumann &
                     .or. bclo(3) .eq. amrex_lo_inflow)) then
                   fac = 2.d0*fac                   
                else if (k .eq. ndhi(3) .and. &
                     (    bchi(3) .eq. amrex_lo_neumann &
                     .or. bchi(3) .eq. amrex_lo_inflow)) then
                   fac = 2.d0*fac
                end if

                resid(i,j,k) = (rhs(i,j,k) - resid(i,j,k)) * fac
             else
                resid(i,j,k) = 0.d0
             end if
          end do
       end do
    end do
  end subroutine amrex_mlndlap_crse_resid


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


  subroutine amrex_mlndlap_zero_fine (lo, hi, phi, dlo, dhi, msk, mlo, mhi, fine_flag) &
       bind(c, name='amrex_mlndlap_zero_fine')
    integer, dimension(3), intent(in) :: lo, hi, dlo, dhi, mlo, mhi
    real(amrex_real), intent(inout) :: phi(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    integer         , intent(in   ) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    integer, intent(in) :: fine_flag

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             ! Testing if the node is covered by a fine level in computing
             ! coarse sync residual
             if (all(msk(i-1:i,j-1:j,k-1:k).eq.fine_flag)) then
                phi(i,j,k) = 0.d0
             end if
          end do
       end do
    end do
  end subroutine amrex_mlndlap_zero_fine


#ifdef AMREX_USE_EB

  subroutine amrex_mlndlap_mknewu_eb (lo, hi, u, ulo, uhi, p, plo, phi, sig, slo, shi, &
       vfrac, vlo, vhi, intg, glo, ghi, dxinv) bind(c,name='amrex_mlndlap_mknewu_eb')
    integer, dimension(3), intent(in) :: lo, hi, ulo, uhi, plo, phi, slo, shi, vlo, vhi, glo, ghi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) ::   u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)
    real(amrex_real), intent(in   ) ::   p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(amrex_real), intent(in   )::vfrac(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
    real(amrex_real), intent(in   ) ::intg(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),n_Sintg)

    integer :: i, j, k
    real(amrex_real) :: dpdx, dpdy, dpdz, dpp_xy, dpp_xz, dpp_yz, dpp_xyz

    do    k = lo(3), hi(3)
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (vfrac(i,j,k) .eq. 0.d0) then
             u(i,j,k,1) = 0.d0
             u(i,j,k,2) = 0.d0
             u(i,j,k,3) = 0.d0
          else
             dpdx = 0.25d0*(-p(i,j,k  )+p(i+1,j,k  )-p(i,j+1,k  )+p(i+1,j+1,k  ) &
                            -p(i,j,k+1)+p(i+1,j,k+1)-p(i,j+1,k+1)+p(i+1,j+1,k+1))
             dpdy = 0.25d0*(-p(i,j,k  )-p(i+1,j,k  )+p(i,j+1,k  )+p(i+1,j+1,k  ) &
                            -p(i,j,k+1)-p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1))
             dpdz = 0.25d0*(-p(i,j,k  )-p(i+1,j,k  )-p(i,j+1,k  )-p(i+1,j+1,k  ) &
                            +p(i,j,k+1)+p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1))

             dpp_xy = (p(i+1,j+1,k+1) - p(i,j+1,k+1) - p(i+1,j,k+1) + p(i,j,k+1) &
                      +p(i+1,j+1,k  ) - p(i,j+1,k  ) - p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k)

             dpp_xz = (p(i+1,j+1,k+1) - p(i,j+1,k+1) + p(i+1,j,k+1) - p(i,j,k+1) &
                      -p(i+1,j+1,k  ) + p(i,j+1,k  ) - p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k)

             dpp_yz = (p(i+1,j+1,k+1) + p(i,j+1,k+1) - p(i+1,j,k+1) - p(i,j,k+1) &
                      -p(i+1,j+1,k  ) - p(i,j+1,k  ) + p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k)

             dpp_xyz = (p(i+1,j+1,k+1) - p(i,j+1,k+1) - p(i+1,j,k+1) + p(i,j,k+1) &
                       -p(i+1,j+1,k  ) + p(i,j+1,k  ) + p(i+1,j,k  ) - p(i,j,k  ) ) / vfrac(i,j,k)

             u(i,j,k,1) = u(i,j,k,1) - sig(i,j,k)*dxinv(1)*(dpdx + .5d0*intg(i,j,k,i_S_y  )*dpp_xy + &
                                                                   .5d0*intg(i,j,k,i_S_z  )*dpp_xz + &
                                                                        intg(i,j,k,i_S_y_z)*dpp_xyz )
             u(i,j,k,2) = u(i,j,k,2) - sig(i,j,k)*dxinv(2)*(dpdy + .5d0*intg(i,j,k,i_S_x  )*dpp_xy + &
                                                                   .5d0*intg(i,j,k,i_S_z  )*dpp_yz + &
                                                                        intg(i,j,k,i_S_x_z)*dpp_xyz )
             u(i,j,k,3) = u(i,j,k,3) - sig(i,j,k)*dxinv(3)*(dpdz + .5d0*intg(i,j,k,i_S_x  )*dpp_xz + &
                                                                   .5d0*intg(i,j,k,i_S_y  )*dpp_yz + & 
                                                                        intg(i,j,k,i_S_x_y)*dpp_xyz )

          end if
       end do
    end do
    end do

  end subroutine amrex_mlndlap_mknewu_eb

  subroutine amrex_mlndlap_rhcc_eb (lo, hi, rhs, rlo, rhi, rhcc, clo, chi, &
       vfrac, flo, fhi, intg, glo, ghi, msk, mlo, mhi) &
       bind(c,name='amrex_mlndlap_rhcc_eb')
    integer, dimension(3) :: lo, hi, rlo, rhi, clo, chi, flo, fhi, glo, ghi, mlo, mhi
    real(amrex_real), intent(inout) :: rhs (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(in   ) :: rhcc(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(amrex_real), intent(in   ) :: vfrac(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(amrex_real), intent(in   ) :: intg(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),n_Sintg)
    integer,          intent(in   ) :: msk (mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .ne. dirichlet) then
                rhs(i,j,k) = &
                     &            rhcc(i  ,j  ,k  ) * &
                     ( 0.125d0 * vfrac(i  ,j  ,k  ) &
                     + 0.25d0 * (-intg(i  ,j  ,k  ,i_S_x) &
                     &           -intg(i  ,j  ,k  ,i_S_y) &
                     &           -intg(i  ,j  ,k  ,i_S_z)) &
                     + 0.5d0 * (  intg(i  ,j  ,k  ,i_S_x_y) &
                     &           +intg(i  ,j  ,k  ,i_S_x_z) &
                     &           +intg(i  ,j  ,k  ,i_S_y_z)) &
                     +         ( -intg(i  ,j  ,k  ,i_S_xyz))) &
                     !
                     +            rhcc(i-1,j  ,k  ) * &
                     ( 0.125d0 * vfrac(i-1,j  ,k  ) &
                     + 0.25d0 * ( intg(i-1,j  ,k  ,i_S_x) &
                     &           -intg(i-1,j  ,k  ,i_S_y) &
                     &           -intg(i-1,j  ,k  ,i_S_z)) &
                     + 0.5d0 * ( -intg(i-1,j  ,k  ,i_S_x_y) &
                     &           -intg(i-1,j  ,k  ,i_S_x_z) &
                     &           +intg(i-1,j  ,k  ,i_S_y_z)) &
                     +         (  intg(i-1,j  ,k  ,i_S_xyz))) &
                     !
                     +            rhcc(i  ,j-1,k  ) * &
                     ( 0.125d0 * vfrac(i  ,j-1,k  ) &
                     + 0.25d0 * (-intg(i  ,j-1,k  ,i_S_x) &
                     &           +intg(i  ,j-1,k  ,i_S_y) &
                     &           -intg(i  ,j-1,k  ,i_S_z)) &
                     + 0.5d0 * ( -intg(i  ,j-1,k  ,i_S_x_y) &
                     &           +intg(i  ,j-1,k  ,i_S_x_z) &
                     &           -intg(i  ,j-1,k  ,i_S_y_z)) &
                     +         (  intg(i  ,j-1,k  ,i_S_xyz))) &
                     !
                     +            rhcc(i-1,j-1,k  ) * &
                     ( 0.125d0 * vfrac(i-1,j-1,k  ) &
                     + 0.25d0 * ( intg(i-1,j-1,k  ,i_S_x) &
                     &           +intg(i-1,j-1,k  ,i_S_y) &
                     &           -intg(i-1,j-1,k  ,i_S_z)) &
                     + 0.5d0 * (  intg(i-1,j-1,k  ,i_S_x_y) &
                     &           -intg(i-1,j-1,k  ,i_S_x_z) &
                     &           -intg(i-1,j-1,k  ,i_S_y_z)) &
                     +         ( -intg(i-1,j-1,k  ,i_S_xyz))) &
                     !
                     +            rhcc(i  ,j  ,k-1) * &
                     ( 0.125d0 * vfrac(i  ,j  ,k-1) &
                     + 0.25d0 * (-intg(i  ,j  ,k-1,i_S_x) &
                     &           -intg(i  ,j  ,k-1,i_S_y) &
                     &           +intg(i  ,j  ,k-1,i_S_z)) &
                     + 0.5d0 * (  intg(i  ,j  ,k-1,i_S_x_y) &
                     &           -intg(i  ,j  ,k-1,i_S_x_z) &
                     &           -intg(i  ,j  ,k-1,i_S_y_z)) &
                     +         (  intg(i  ,j  ,k-1,i_S_xyz))) &
                     !
                     +            rhcc(i-1,j  ,k-1) * &
                     ( 0.125d0 * vfrac(i-1,j  ,k-1) &
                     + 0.25d0 * ( intg(i-1,j  ,k-1,i_S_x) &
                     &           -intg(i-1,j  ,k-1,i_S_y) &
                     &           +intg(i-1,j  ,k-1,i_S_z)) &
                     + 0.5d0 * ( -intg(i-1,j  ,k-1,i_S_x_y) &
                     &           +intg(i-1,j  ,k-1,i_S_x_z) &
                     &           -intg(i-1,j  ,k-1,i_S_y_z)) &
                     +         ( -intg(i-1,j  ,k-1,i_S_xyz))) &
                     !
                     +            rhcc(i  ,j-1,k-1) * &
                     ( 0.125d0 * vfrac(i  ,j-1,k-1) &
                     + 0.25d0 * (-intg(i  ,j-1,k-1,i_S_x) &
                     &           +intg(i  ,j-1,k-1,i_S_y) &
                     &           +intg(i  ,j-1,k-1,i_S_z)) &
                     + 0.5d0 * ( -intg(i  ,j-1,k-1,i_S_x_y) &
                     &           -intg(i  ,j-1,k-1,i_S_x_z) &
                     &           +intg(i  ,j-1,k-1,i_S_y_z)) &
                     +         ( -intg(i  ,j-1,k-1,i_S_xyz))) &
                     !
                     +            rhcc(i-1,j-1,k-1) * &
                     ( 0.125d0 * vfrac(i-1,j-1,k-1) &
                     + 0.25d0 * ( intg(i-1,j-1,k-1,i_S_x) &
                     &           +intg(i-1,j-1,k-1,i_S_y) &
                     &           +intg(i-1,j-1,k-1,i_S_z)) &
                     + 0.5d0 * (  intg(i-1,j-1,k-1,i_S_x_y) &
                     &           +intg(i-1,j-1,k-1,i_S_x_z) &
                     &           +intg(i-1,j-1,k-1,i_S_y_z)) &
                     +         (  intg(i-1,j-1,k-1,i_S_xyz)))
             else
                rhs(i,j,k) = 0.d0
             end if
          end do
       end do
    end do
  end subroutine amrex_mlndlap_rhcc_eb

#endif

end module amrex_mlnodelap_3d_module
