module amrex_mlnodelap_2d_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  use amrex_lo_bctypes_module, only : amrex_lo_dirichlet, amrex_lo_neumann, amrex_lo_inflow, amrex_lo_periodic
  implicit none

  ! external dirichlet at physical boundary or internal dirichlet at crse/fine boundary
  integer, parameter :: dirichlet = 1

  private
  public :: amrex_mlndlap_avgdown_coeff, amrex_mlndlap_fillbc_cc, amrex_mlndlap_divu, &
       amrex_mlndlap_applybc, amrex_mlndlap_restriction, amrex_mlndlap_mknewu, &
       amrex_mlndlap_adotx_ha, amrex_mlndlap_adotx_aa, &
       amrex_mlndlap_jacobi_ha, amrex_mlndlap_jacobi_aa, &
       amrex_mlndlap_gauss_seidel_ha, amrex_mlndlap_gauss_seidel_aa, &
       amrex_mlndlap_interpolation_ha, amrex_mlndlap_interpolation_aa, &
       amrex_mlndlap_zero_fine, amrex_mlndlap_crse_resid, amrex_mlndlap_any_zero, &
       amrex_mlndlap_set_dirichlet_mask, amrex_mlndlap_divu_fine_contrib, &
       amrex_mlndlap_divu_cf_contrib, amrex_mlndlap_res_fine_contrib, &
       amrex_mlndlap_res_cf_contrib, amrex_mlndlap_fixup_res_mask

contains

  subroutine amrex_mlndlap_avgdown_coeff (lo, hi, crse, clo, chi, fine, flo, fhi, idim) &
       bind(c,name='amrex_mlndlap_avgdown_coeff')
    integer, dimension(2), intent(in) :: lo, hi, clo, chi, flo, fhi
    integer, intent(in) :: idim
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2))
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1),flo(2):fhi(2))

    integer :: i,j

    if (idim .eq. 0) then
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             crse(i,j) = (fine(2*i,2*j)+fine(2*i,2*j+1))*(fine(2*i+1,2*j)+fine(2*i+1,2*j+1))/ &
                  (fine(2*i,2*j)+fine(2*i+1,2*j)+fine(2*i,2*j+1)+fine(2*i+1,2*j+1))
          end do
       end do
    else
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             crse(i,j) = (fine(2*i,2*j)+fine(2*i+1,2*j))*(fine(2*i,2*j+1)+fine(2*i+1,2*j+1))/ &
                  (fine(2*i,2*j)+fine(2*i+1,2*j)+fine(2*i,2*j+1)+fine(2*i+1,2*j+1))
          end do
       end do
    end if

  end subroutine amrex_mlndlap_avgdown_coeff


  subroutine amrex_mlndlap_fillbc_cc (sigma, slo, shi, dlo, dhi, bclo, bchi) &
       bind(c, name='amrex_mlndlap_fillbc_cc')
    integer, dimension(2), intent(in) :: slo, shi, dlo, dhi, bclo, bchi
    real(amrex_real), intent(inout) :: sigma(slo(1):shi(1),slo(2):shi(2))

    integer :: ilo, ihi, jlo, jhi

    ilo = max(dlo(1), slo(1))
    ihi = min(dhi(1), shi(1))
    jlo = max(dlo(2), slo(2))
    jhi = min(dhi(2), shi(2))

    if (bclo(1) .ne. amrex_lo_periodic .and. slo(1) .lt. dlo(1)) then
       sigma(dlo(1)-1,jlo:jhi) = sigma(dlo(1),jlo:jhi)
    end if
    
    if (bchi(1) .ne. amrex_lo_periodic .and. shi(1) .gt. dhi(1)) then
       sigma(dhi(1)+1,jlo:jhi) = sigma(dhi(1),jlo:jhi)
    end if

    if (bclo(2) .ne. amrex_lo_periodic .and. slo(2) .lt. dlo(2)) then
       sigma(ilo:ihi,dlo(2)-1) = sigma(ilo:ihi,dlo(2))
    end if

    if (bchi(2) .ne. amrex_lo_periodic .and. shi(2) .gt. dhi(2)) then
       sigma(ilo:ihi,dhi(2)+1) = sigma(ilo:ihi,dhi(2))
    end if

    if (slo(1) .lt. dlo(1) .and. slo(2) .lt. dlo(2)) then
       if (bclo(1) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dlo(2)-1) = sigma(dlo(1),dlo(2)-1)
       else if (bclo(2) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dlo(2)-1) = sigma(dlo(1)-1,dlo(2))
       end if
    end if

    if (shi(1) .gt. dhi(1) .and. slo(2) .lt. dlo(2)) then
       if (bchi(1) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dlo(2)-1) = sigma(dhi(1),dlo(2)-1)
       else if (bclo(2) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dlo(2)-1) = sigma(dhi(1)+1,dlo(2))
       end if
    end if

    if (slo(1) .lt. dlo(1) .and. shi(2) .gt. dhi(2)) then
       if (bclo(1) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dhi(2)+1) = sigma(dlo(1),dhi(2)+1)
       else if (bchi(2) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dhi(2)+1) = sigma(dlo(1)-1,dhi(2))
       end if
    end if

    if (shi(1) .gt. dhi(1) .and. shi(2) .gt. dhi(2)) then
       if (bchi(1) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dhi(2)+1) = sigma(dhi(1),dhi(2)+1)
       else if (bchi(2) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dhi(2)+1) = sigma(dhi(1)+1,dhi(2))
       end if
    end if

  end subroutine amrex_mlndlap_fillbc_cc


  subroutine amrex_mlndlap_divu (lo, hi, rhs, rlo, rhi, vel, vlo, vhi, msk, mlo, mhi, &
       dxinv, ndlo, ndhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_divu')
    integer, dimension(2), intent(in) :: lo, hi, rlo, rhi, vlo, vhi, mlo, mhi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1),vlo(2):vhi(2),2)
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))

    integer :: i,j
    real(amrex_real) :: facx, facy

    facx = 0.5d0*dxinv(1)
    facy = 0.5d0*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (msk(i,j) .ne. dirichlet) then
             rhs(i,j) = facx*(-vel(i-1,j-1,1)+vel(i,j-1,1)-vel(i-1,j,1)+vel(i,j,1)) &
                  &   + facy*(-vel(i-1,j-1,2)-vel(i,j-1,2)+vel(i-1,j,2)+vel(i,j,2))
          else
             rhs(i,j) = 0.d0
          end if
       end do
    end do

    if (lo(1) .eq. ndlo(1)) then
       if (bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) then 
          rhs(lo(1),lo(2):hi(2)) = 2.d0*rhs(lo(1),lo(2):hi(2))
       end if
    end if

    if (hi(1) .eq. ndhi(1)) then
       if (bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) then
          rhs(hi(1),lo(2):hi(2)) = 2.d0*rhs(hi(1),lo(2):hi(2))
       end if
    end if

    if (lo(2) .eq. ndlo(2)) then
       if (bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) then
          rhs(lo(1):hi(1),lo(2)) = 2.d0*rhs(lo(1):hi(1),lo(2))
       end if
    end if

    if (hi(2) .eq. ndhi(2)) then
       if (bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) then
          rhs(lo(1):hi(1),hi(2)) = 2.d0*rhs(lo(1):hi(1),hi(2))
       end if
    end if

  end subroutine amrex_mlndlap_divu


  subroutine amrex_mlndlap_add_rhcc (lo, hi, rhs, rlo, rhi, rhcc, clo, chi, msk, mlo, mhi) &
       bind(c,name='amrex_mlndlap_add_rhcc')
    integer, dimension(2) :: lo, hi, rlo, rhi, clo, chi, mlo, mhi
    real(amrex_real), intent(inout) :: rhs (rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(in   ) :: rhcc(clo(1):chi(1),clo(2):chi(2))
    integer,          intent(in   ) :: msk (mlo(1):mhi(1),mlo(2):mhi(2))

    integer :: i,j

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (msk(i,j) .ne. dirichlet) then
             rhs(i,j) = rhs(i,j) + 0.25d0*(rhcc(i-1,j-1)+rhcc(i,j-1)+rhcc(i-1,j)+rhcc(i,j))
          end if
       end do
    end do
  end subroutine amrex_mlndlap_add_rhcc


  subroutine amrex_mlndlap_mknewu (lo, hi, u, ulo, uhi, p, plo, phi, sig, slo, shi, dxinv) &
       bind(c,name='amrex_mlndlap_mknewu')
    integer, dimension(2), intent(in) :: lo, hi, ulo, uhi, plo, phi, slo, shi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) ::   u(ulo(1):uhi(1),ulo(2):uhi(2),2)
    real(amrex_real), intent(in   ) ::   p(plo(1):phi(1),plo(2):phi(2))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1),slo(2):shi(2))

    integer :: i, j
    real(amrex_real) :: facx, facy

    facx = 0.5d0*dxinv(1)
    facy = 0.5d0*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          u(i,j,1) = u(i,j,1) - sig(i,j)*facx*(-p(i,j)+p(i+1,j)-p(i,j+1)+p(i+1,j+1))
          u(i,j,2) = u(i,j,2) - sig(i,j)*facy*(-p(i,j)-p(i+1,j)+p(i,j+1)+p(i+1,j+1))
       end do
    end do
  end subroutine amrex_mlndlap_mknewu


  subroutine amrex_mlndlap_applybc (phi, hlo, hhi, dlo, dhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_applybc')
    integer, dimension(2) :: hlo, hhi, dlo, dhi, bclo, bchi
    real(amrex_real), intent(inout) :: phi(hlo(1):hhi(1),hlo(2):hhi(2))

    integer :: ilo, ihi, jlo, jhi

    if (any(bclo.eq.amrex_lo_inflow) .or. any(bchi.eq.amrex_lo_inflow)) then
       call amrex_error("amrex_mlndlap_applybc: inflow not supported yet")
    end if

    ilo = max(dlo(1), hlo(1))
    ihi = min(dhi(1), hhi(1))
    jlo = max(dlo(2), hlo(2))
    jhi = min(dhi(2), hhi(2))

    ! neumann

    if (bclo(1) .eq. amrex_lo_neumann .and. hlo(1) .lt. dlo(1)) then
       phi(dlo(1)-1,jlo:jhi) = phi(dlo(1)+1,jlo:jhi)
    end if

    if (bchi(1) .eq. amrex_lo_neumann .and. hhi(1) .gt. dhi(1)) then
       phi(dhi(1)+1,jlo:jhi) = phi(dhi(1)-1,jlo:jhi)
    end if

    if (bclo(2) .eq. amrex_lo_neumann .and. hlo(2) .lt. dlo(2)) then
       phi(ilo:ihi,dlo(2)-1) = phi(ilo:ihi,dlo(2)+1)
    end if

    if (bchi(2) .eq. amrex_lo_neumann .and. hhi(2) .gt. dhi(2)) then
       phi(ilo:ihi,dhi(2)+1) = phi(ilo:ihi,dhi(2)-1)
    end if

    ! corners

    if (hlo(1) .lt. dlo(1) .and. hlo(2) .lt. dlo(2)) then
       if (bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dlo(2)-1) = phi(dlo(1)+1,dlo(2)-1)
       else if (bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dlo(2)-1) = phi(dlo(1)-1,dlo(2)+1)
       end if
    end if

    if (hhi(1) .gt. dhi(1) .and. hlo(2) .lt. dlo(2)) then
       if (bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dlo(2)-1) = phi(dhi(1)-1,dlo(2)-1)
       else if (bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dlo(2)-1) = phi(dhi(1)+1,dlo(2)+1)
       end if
    end if

    if (hlo(1) .lt. dlo(1) .and. hhi(2) .gt. dhi(2)) then
       if (bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dhi(2)+1) = phi(dlo(1)+1,dhi(2)+1)
       else if (bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dhi(2)+1) = phi(dlo(1)-1,dhi(2)-1)
       end if
    end if

    if (hhi(1) .gt. dhi(1) .and. hhi(2) .gt. dhi(2)) then
       if (bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dhi(2)+1) = phi(dhi(1)-1,dhi(2)+1)
       else  if (bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dhi(2)+1) = phi(dhi(1)+1,dhi(2)-1)
       end if
    end if

  end subroutine amrex_mlndlap_applybc


  subroutine amrex_mlndlap_adotx_ha (lo, hi, y, ylo, yhi, x, xlo, xhi, &
       sx, sxlo, sxhi, sy, sylo, syhi, dg, dlo, dhi, msk, mlo, mhi, dxinv, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_adotx_ha')
    integer, dimension(2), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, sxlo, sxhi, sylo, syhi, dlo, dhi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) ::  y( ylo(1): yhi(1), ylo(2): yhi(2))
    real(amrex_real), intent(in   ) ::  x( xlo(1): xhi(1), xlo(2): xhi(2))
    real(amrex_real), intent(in   ) :: sx(sxlo(1):sxhi(1),sxlo(2):sxhi(2))
    real(amrex_real), intent(in   ) :: sy(sylo(1):syhi(1),sylo(2):syhi(2))
    real(amrex_real)                :: dg( dlo(1): dhi(1), dlo(2): dhi(2), 2)
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))
 
    integer :: i,j
    real(amrex_real) :: facx, facy

    facx = (1.d0/6.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/6.d0)*dxinv(2)*dxinv(2)

    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)-1, hi(1)
          dg(i,j,1) = x(i+1,j) - x(i,j)
       end do
    end do

    do    j = lo(2)-1, hi(2)
       do i = lo(1)-1, hi(1)+1
          dg(i,j,2) = x(i,j+1) - x(i,j)
       end do
    end do

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (msk(i,j) .ne. dirichlet) then
             y(i,j) = facx*(-sx(i-1,j-1)*(dg(i-1,j-1,1)+2.d0*dg(i-1,j  ,1)) &
                  &         +sx(i  ,j-1)*(dg(i  ,j-1,1)+2.d0*dg(i  ,j  ,1)) &
                  &         -sx(i-1,j  )*(dg(i-1,j+1,1)+2.d0*dg(i-1,j  ,1)) &
                  &         +sx(i  ,j  )*(dg(i  ,j+1,1)+2.d0*dg(i  ,j  ,1))) &
                  +   facy*(-sy(i-1,j-1)*(dg(i-1,j-1,2)+2.d0*dg(i  ,j-1,2)) &
                  &         -sy(i  ,j-1)*(dg(i+1,j-1,2)+2.d0*dg(i  ,j-1,2)) &
                  &         +sy(i-1,j  )*(dg(i-1,j  ,2)+2.d0*dg(i  ,j  ,2)) &
                  &         +sy(i  ,j  )*(dg(i+1,j  ,2)+2.d0*dg(i  ,j  ,2)))
          else
             y(i,j) = 0.d0
          end if
       end do
    end do

  end subroutine amrex_mlndlap_adotx_ha


  subroutine amrex_mlndlap_adotx_aa (lo, hi, y, ylo, yhi, x, xlo, xhi, &
       sig, slo, shi, dg, dlo, dhi, msk, mlo, mhi, dxinv, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_adotx_aa')
    integer, dimension(2), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, slo, shi, dlo, dhi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) ::   y(ylo(1):yhi(1),ylo(2):yhi(2))
    real(amrex_real), intent(in   ) ::   x(xlo(1):xhi(1),xlo(2):xhi(2))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1),slo(2):shi(2))
    real(amrex_real)                ::  dg(dlo(1):dhi(1),dlo(2):dhi(2), 2)
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))

    integer :: i,j
    real(amrex_real) :: facx, facy

    facx = (1.d0/6.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/6.d0)*dxinv(2)*dxinv(2)

    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)-1, hi(1)
          dg(i,j,1) = x(i+1,j) - x(i,j)
       end do
    end do

    do    j = lo(2)-1, hi(2)
       do i = lo(1)-1, hi(1)+1
          dg(i,j,2) = x(i,j+1) - x(i,j)
       end do
    end do

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (msk(i,j) .ne. dirichlet) then
             y(i,j) = facx*(-sig(i-1,j-1)*(dg(i-1,j-1,1)+2.d0*dg(i-1,j  ,1)) &
                  &         +sig(i  ,j-1)*(dg(i  ,j-1,1)+2.d0*dg(i  ,j  ,1)) &
                  &         -sig(i-1,j  )*(dg(i-1,j+1,1)+2.d0*dg(i-1,j  ,1)) &
                  &         +sig(i  ,j  )*(dg(i  ,j+1,1)+2.d0*dg(i  ,j  ,1))) &
                  +   facy*(-sig(i-1,j-1)*(dg(i-1,j-1,2)+2.d0*dg(i  ,j-1,2)) &
                  &         -sig(i  ,j-1)*(dg(i+1,j-1,2)+2.d0*dg(i  ,j-1,2)) &
                  &         +sig(i-1,j  )*(dg(i-1,j  ,2)+2.d0*dg(i  ,j  ,2)) &
                  &         +sig(i  ,j  )*(dg(i+1,j  ,2)+2.d0*dg(i  ,j  ,2)))
          else
             y(i,j) = 0.d0
          end if
       end do
    end do

  end subroutine amrex_mlndlap_adotx_aa


  subroutine amrex_mlndlap_jacobi_ha (lo, hi, sol, slo, shi, Ax, alo, ahi, rhs, rlo, rhi, &
       sx, sxlo, sxhi, sy, sylo, syhi, msk, mlo, mhi, dxinv, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_jacobi_ha')
    integer, dimension(2),intent(in) :: lo,hi,slo,shi,alo,ahi,rlo,rhi,sxlo,sxhi,sylo,syhi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: sol( slo(1): shi(1), slo(2): shi(2))
    real(amrex_real), intent(in   ) :: Ax ( alo(1): ahi(1), alo(2): ahi(2))
    real(amrex_real), intent(in   ) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2))
    real(amrex_real), intent(in   ) :: sx (sxlo(1):sxhi(1),sxlo(2):sxhi(2))
    real(amrex_real), intent(in   ) :: sy (sylo(1):syhi(1),sylo(2):syhi(2))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))

    integer :: i,j
    real(amrex_real) :: facx, facy
    real(amrex_real), parameter :: omega = 2.d0/3.d0

    facx = -2.d0 * (1.d0/6.d0)*dxinv(1)*dxinv(1)
    facy = -2.d0 * (1.d0/6.d0)*dxinv(2)*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (msk(i,j) .ne. dirichlet) then
             sol(i,j) = sol(i,j) + omega * (rhs(i,j) - Ax(i,j)) &
                  / (facx*(sx(i-1,j-1)+sx(i,j-1)+sx(i-1,j)+sx(i,j)) &
                  +  facy*(sy(i-1,j-1)+sy(i,j-1)+sy(i-1,j)+sy(i,j)))
          else
             sol(i,j) = 0.d0
          end if
       end do
    end do

  end subroutine amrex_mlndlap_jacobi_ha


  subroutine amrex_mlndlap_jacobi_aa (lo, hi, sol, slo, shi, Ax, alo, ahi, rhs, rlo, rhi, &
       sig, sglo, sghi, msk, mlo, mhi, dxinv, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_jacobi_aa')
    integer, dimension(2),intent(in) :: lo,hi,slo,shi,alo,ahi,rlo,rhi,sglo,sghi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: sol( slo(1): shi(1), slo(2): shi(2))
    real(amrex_real), intent(in   ) :: Ax ( alo(1): ahi(1), alo(2): ahi(2))
    real(amrex_real), intent(in   ) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2))
    real(amrex_real), intent(in   ) :: sig(sglo(1):sghi(1),sglo(2):sghi(2))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))

    integer :: i,j
    real(amrex_real) :: facx, facy, fac
    real(amrex_real), parameter :: omega = 2.d0/3.d0

    facx = -2.d0 * (1.d0/6.d0)*dxinv(1)*dxinv(1)
    facy = -2.d0 * (1.d0/6.d0)*dxinv(2)*dxinv(2)
    fac = facx + facy

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (msk(i,j) .ne. dirichlet) then
             sol(i,j) = sol(i,j) + omega * (rhs(i,j) - Ax(i,j)) &
                  / (fac*(sig(i-1,j-1)+sig(i,j-1)+sig(i-1,j)+sig(i,j)))
          else
             sol(i,j) = 0.d0
          end if
       end do
    end do

  end subroutine amrex_mlndlap_jacobi_aa


  subroutine amrex_mlndlap_gauss_seidel_ha (lo, hi, sol, slo, shi, rhs, rlo, rhi, &
       sx, sxlo, sxhi, sy, sylo, syhi, msk, mlo, mhi, dxinv, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_gauss_seidel_ha')
    integer, dimension(2),intent(in) :: lo,hi,slo,shi,rlo,rhi,sxlo,sxhi,sylo,syhi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: sol( slo(1): shi(1), slo(2): shi(2))
    real(amrex_real), intent(in   ) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2))
    real(amrex_real), intent(in   ) :: sx (sxlo(1):sxhi(1),sxlo(2):sxhi(2))
    real(amrex_real), intent(in   ) :: sy (sylo(1):syhi(1),sylo(2):syhi(2))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))

    integer :: i,j
    real(amrex_real) :: facx, facy, Ax
    real(amrex_real) :: dgx_mm, dgx_m0, dgx_mp, dgx_0m, dgx_00, dgx_0p
    real(amrex_real) :: dgy_mm, dgy_m0, dgy_0m, dgy_00, dgy_pm, dgy_p0

    facx = (1.d0/6.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/6.d0)*dxinv(2)*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (msk(i,j) .ne. dirichlet) then
             dgx_mm = sol(i  ,j-1) - sol(i-1,j-1)
             dgx_m0 = sol(i  ,j  ) - sol(i-1,j  )
             dgx_mp = sol(i  ,j+1) - sol(i-1,j+1)
             dgx_0m = sol(i+1,j-1) - sol(i  ,j-1)
             dgx_00 = sol(i+1,j  ) - sol(i  ,j  )
             dgx_0p = sol(i+1,j+1) - sol(i  ,j+1)

             dgy_mm = sol(i-1,j  ) - sol(i-1,j-1)
             dgy_m0 = sol(i-1,j+1) - sol(i-1,j  )
             dgy_0m = sol(i  ,j  ) - sol(i  ,j-1)
             dgy_00 = sol(i  ,j+1) - sol(i  ,j  )
             dgy_pm = sol(i+1,j  ) - sol(i+1,j-1)
             dgy_p0 = sol(i+1,j+1) - sol(i+1,j  )

             Ax =   facx*(-sx(i-1,j-1)*(dgx_mm + 2.d0*dgx_m0) &
                  &       +sx(i  ,j-1)*(dgx_0m + 2.d0*dgx_00) &
                  &       -sx(i-1,j  )*(dgx_mp + 2.d0*dgx_m0) &
                  &       +sx(i  ,j  )*(dgx_0p + 2.d0*dgx_00)) &
                  + facy*(-sy(i-1,j-1)*(dgy_mm + 2.d0*dgy_0m) &
                  &       -sy(i  ,j-1)*(dgy_pm + 2.d0*dgy_0m) &
                  &       +sy(i-1,j  )*(dgy_m0 + 2.d0*dgy_00) &
                  &       +sy(i  ,j  )*(dgy_p0 + 2.d0*dgy_00))
             
             sol(i,j) = sol(i,j) + (rhs(i,j) - Ax) * (-0.5d0) &
                  / (facx*(sx(i-1,j-1)+sx(i,j-1)+sx(i-1,j)+sx(i,j)) &
                  +  facy*(sy(i-1,j-1)+sy(i,j-1)+sy(i-1,j)+sy(i,j)))
          else
             sol(i,j) = 0.d0
          end if
       end do
    end do

  end subroutine amrex_mlndlap_gauss_seidel_ha


  subroutine amrex_mlndlap_gauss_seidel_aa (lo, hi, sol, slo, shi, rhs, rlo, rhi, &
       sig, sglo, sghi, msk, mlo, mhi, dxinv, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_gauss_seidel_aa')
    integer, dimension(2),intent(in) :: lo,hi,slo,shi,rlo,rhi,sglo,sghi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: sol( slo(1): shi(1), slo(2): shi(2))
    real(amrex_real), intent(in   ) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2))
    real(amrex_real), intent(in   ) :: sig(sglo(1):sghi(1),sglo(2):sghi(2))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))

    integer :: i,j
    real(amrex_real) :: facx, facy, Ax, fac
    real(amrex_real) :: dgx_mm, dgx_m0, dgx_mp, dgx_0m, dgx_00, dgx_0p
    real(amrex_real) :: dgy_mm, dgy_m0, dgy_0m, dgy_00, dgy_pm, dgy_p0

    facx = (1.d0/6.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/6.d0)*dxinv(2)*dxinv(2)
    fac = facx + facy

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (msk(i,j) .ne. dirichlet) then
             dgx_mm = sol(i  ,j-1) - sol(i-1,j-1)
             dgx_m0 = sol(i  ,j  ) - sol(i-1,j  )
             dgx_mp = sol(i  ,j+1) - sol(i-1,j+1)
             dgx_0m = sol(i+1,j-1) - sol(i  ,j-1)
             dgx_00 = sol(i+1,j  ) - sol(i  ,j  )
             dgx_0p = sol(i+1,j+1) - sol(i  ,j+1)
             
             dgy_mm = sol(i-1,j  ) - sol(i-1,j-1)
             dgy_m0 = sol(i-1,j+1) - sol(i-1,j  )
             dgy_0m = sol(i  ,j  ) - sol(i  ,j-1)
             dgy_00 = sol(i  ,j+1) - sol(i  ,j  )
             dgy_pm = sol(i+1,j  ) - sol(i+1,j-1)
             dgy_p0 = sol(i+1,j+1) - sol(i+1,j  )
             
             Ax =   facx*(-sig(i-1,j-1)*(dgx_mm + 2.d0*dgx_m0) &
                  &       +sig(i  ,j-1)*(dgx_0m + 2.d0*dgx_00) &
                  &       -sig(i-1,j  )*(dgx_mp + 2.d0*dgx_m0) &
                  &       +sig(i  ,j  )*(dgx_0p + 2.d0*dgx_00)) &
                  + facy*(-sig(i-1,j-1)*(dgy_mm + 2.d0*dgy_0m) &
                  &       -sig(i  ,j-1)*(dgy_pm + 2.d0*dgy_0m) &
                  &       +sig(i-1,j  )*(dgy_m0 + 2.d0*dgy_00) &
                  &       +sig(i  ,j  )*(dgy_p0 + 2.d0*dgy_00))
             
             sol(i,j) = sol(i,j) + (rhs(i,j) - Ax) * (-0.5d0) &
                  / (fac*(sig(i-1,j-1)+sig(i,j-1)+sig(i-1,j)+sig(i,j)))
          else
             sol(i,j) = 0.d0
          end if
       end do
    end do

  end subroutine amrex_mlndlap_gauss_seidel_aa


  subroutine amrex_mlndlap_restriction (lo, hi, crse, clo, chi, fine, flo, fhi, msk, mlo, mhi, &
       domlo, domhi, bclo, bchi) bind(c,name='amrex_mlndlap_restriction')
    integer, dimension(2), intent(in) :: lo, hi, clo, chi, flo, fhi, mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2))
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1),flo(2):fhi(2))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))

    integer :: i,j, ii, jj
    real(amrex_real), parameter :: fac = 1.d0/16.d0

    do    j = lo(2), hi(2)
       jj = 2*j
       do i = lo(1), hi(1)
          ii = 2*i
          if (msk(ii,jj) .ne. dirichlet) then
             crse(i,j) = fac*(fine(ii-1,jj-1) + 2.d0*fine(ii  ,jj-1) +      fine(ii+1,jj-1) &
                  +      2.d0*fine(ii-1,jj  ) + 4.d0*fine(ii  ,jj  ) + 2.d0*fine(ii+1,jj  ) &
                  +           fine(ii-1,jj+1) + 2.d0*fine(ii  ,jj+1) +      fine(ii+1,jj+1))
          else
             crse(i,j) = 0.d0
          end if
       end do
    end do
  end subroutine amrex_mlndlap_restriction


  subroutine amrex_mlndlap_interpolation_ha (clo, chi, fine, fflo, ffhi, crse, cflo, cfhi, &
       sigx, sxlo, sxhi, sigy, sylo, syhi, msk, mlo, mhi, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_interpolation_ha')
    integer, dimension(2), intent(in) :: clo,chi,fflo,ffhi,cflo,cfhi,sxlo,sxhi,sylo,syhi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in   ) :: crse(cflo(1):cfhi(1),cflo(2):cfhi(2))
    real(amrex_real), intent(inout) :: fine(fflo(1):ffhi(1),fflo(2):ffhi(2))
    real(amrex_real), intent(in   ) :: sigx(sxlo(1):sxhi(1),sxlo(2):sxhi(2))
    real(amrex_real), intent(in   ) :: sigy(sylo(1):syhi(1),sylo(2):syhi(2))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))

    integer :: flo(2), fhi(2), i,j
    logical :: interpx
    real(amrex_real) :: wxm, wxp, wym, wyp

    flo = 2*clo
    fhi = 2*chi

    do    j = clo(2), chi(2)
       do i = clo(1), chi(1)
          if (msk(2*i,2*j) .ne. dirichlet) then
             fine(2*i,2*j) = crse(i,j)
          else
             fine(2*i,2*j) = 0.d0
          end if
       end do
    end do

    interpx = .false.
    do j = flo(2), fhi(2)
       interpx = .not.interpx
       if (interpx) then
          ! interp in x-direction
          do i = flo(1)+1, fhi(1), 2
             if (msk(i,j) .ne. dirichlet) then
                wxm = sigx(i-1,j-1) + sigx(i-1,j)
                wxp = sigx(i  ,j-1) + sigx(i  ,j)
                fine(i,j) = (wxm*fine(i-1,j) + wxp*fine(i+1,j)) / (wxm+wxp)
             else
                fine(i,j) = 0.d0
             end if
          end do
       else
          ! interp in y-direction
          do i = flo(1), fhi(1), 2
             if (msk(i,j) .ne. dirichlet) then
                wym = sigy(i-1,j-1) + sigy(i,j-1)
                wyp = sigy(i-1,j  ) + sigy(i,j  )
                fine(i,j) = (wym*fine(i,j-1) + wyp*fine(i,j+1)) / (wym+wyp)
             else
                fine(i,j) = 0.d0
             end if
          end do
       end if
    end do

    do    j = flo(2)+1, fhi(2), 2
       do i = flo(1)+1, fhi(1), 2
          wxm = sigx(i-1,j-1) + sigx(i-1,j  )
          wxp = sigx(i  ,j-1) + sigx(i  ,j  )
          wym = sigy(i-1,j-1) + sigy(i  ,j-1)
          wyp = sigy(i-1,j  ) + sigy(i  ,j  )
          fine(i,j) = (wxm*fine(i-1,j) + wxp*fine(i+1,j) + wym*fine(i,j-1) + wyp*fine(i,j+1)) &
               / (wxm+wxp+wym+wyp)
       end do
    end do

  end subroutine amrex_mlndlap_interpolation_ha


  subroutine amrex_mlndlap_interpolation_aa (clo, chi, fine, fflo, ffhi, crse, cflo, cfhi, &
       sig, sglo, sghi, msk, mlo, mhi, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_interpolation_aa')
    integer, dimension(2), intent(in) :: clo,chi,fflo,ffhi,cflo,cfhi,sglo,sghi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in   ) :: crse(cflo(1):cfhi(1),cflo(2):cfhi(2))
    real(amrex_real), intent(inout) :: fine(fflo(1):ffhi(1),fflo(2):ffhi(2))
    real(amrex_real), intent(in   ) :: sig (sglo(1):sghi(1),sglo(2):sghi(2))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))

    integer :: flo(2), fhi(2), i,j
    logical :: interpx
    real(amrex_real) :: wxm, wxp, wym, wyp

    flo = 2*clo
    fhi = 2*chi

    do    j = clo(2), chi(2)
       do i = clo(1), chi(1)
          if (msk(2*i,2*j) .ne. dirichlet) then
             fine(2*i,2*j) = crse(i,j)
          else
             fine(2*i,2*j) = 0.d0
          end if
       end do
    end do

    interpx = .false.
    do j = flo(2), fhi(2)
       interpx = .not.interpx
       if (interpx) then
          ! interp in x-direction
          do i = flo(1)+1, fhi(1), 2
             if (msk(i,j) .ne. dirichlet) then
                wxm = sig(i-1,j-1) + sig(i-1,j)
                wxp = sig(i  ,j-1) + sig(i  ,j)
                fine(i,j) = (wxm*fine(i-1,j) + wxp*fine(i+1,j)) / (wxm+wxp)
             else
                fine(i,j) = 0.d0
             end if
          end do
       else
          ! interp in y-direction
          do i = flo(1), fhi(1), 2
             if (msk(i,j) .ne. dirichlet) then
                wym = sig(i-1,j-1) + sig(i,j-1)
                wyp = sig(i-1,j  ) + sig(i,j  )
                fine(i,j) = (wym*fine(i,j-1) + wyp*fine(i,j+1)) / (wym+wyp)
             else
                fine(i,j) = 0.d0
             end if
          end do
       end if
    end do

    do    j = flo(2)+1, fhi(2), 2
       do i = flo(1)+1, fhi(1), 2
          wxm = sig(i-1,j-1) + sig(i-1,j  )
          wxp = sig(i  ,j-1) + sig(i  ,j  )
          wym = sig(i-1,j-1) + sig(i  ,j-1)
          wyp = sig(i-1,j  ) + sig(i  ,j  )
          fine(i,j) = (wxm*fine(i-1,j) + wxp*fine(i+1,j) + wym*fine(i,j-1) + wyp*fine(i,j+1)) &
               / (wxm+wxp+wym+wyp)
       end do
    end do

  end subroutine amrex_mlndlap_interpolation_aa


  subroutine amrex_mlndlap_zero_fine (lo, hi, phi, philo, phihi, msk, mlo, mhi) &
       bind(c, name='amrex_mlndlap_zero_fine')
    integer, dimension(2), intent(in) :: lo, hi, philo, phihi, mlo, mhi
    real(amrex_real), intent(inout) :: phi(philo(1):phihi(1),philo(2):phihi(2))
    integer         , intent(in   ) :: msk(  mlo(1):  mhi(1),  mlo(2):  mhi(2))

    integer :: i,j

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (all(msk(i-1:i,j-1:j).eq.0)) then
             phi(i,j) = 0.d0
          end if
       end do
    end do
  end subroutine amrex_mlndlap_zero_fine


  subroutine amrex_mlndlap_crse_resid (lo, hi, resid, rslo, rshi, rhs, rhlo, rhhi, msk, mlo, mhi) &
       bind(c, name='amrex_mlndlap_crse_resid')
    integer, dimension(2), intent(in) :: lo, hi, rslo, rshi, rhlo, rhhi, mlo, mhi
    real(amrex_real), intent(inout) :: resid(rslo(1):rshi(1),rslo(2):rshi(2))
    real(amrex_real), intent(in   ) :: rhs  (rhlo(1):rhhi(1),rhlo(2):rhhi(2))
    integer         , intent(in   ) :: msk  ( mlo(1): mhi(1), mlo(2): mhi(2))

    integer :: i,j

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (any(msk(i-1:i,j-1:j).eq.0) .and. any(msk(i-1:i,j-1:j).eq.1)) then
             resid(i,j) = rhs(i,j) - resid(i,j)
          else
             resid(i,j) = 0.d0
          end if
       end do
    end do
  end subroutine amrex_mlndlap_crse_resid

  
  function amrex_mlndlap_any_zero (lo, hi, msk, mlo, mhi) result(r) &
       bind(c,name='amrex_mlndlap_any_zero')
    integer, dimension(2), intent(in) :: lo, hi, mlo, mhi
    integer, intent(in   ) :: msk  ( mlo(1): mhi(1), mlo(2): mhi(2))
    integer :: r

    integer :: i,j

    r = 0

    do j = lo(2), hi(2)
       if (r.eq.1) exit
       do i = lo(1), hi(1)
          if (r.eq.1) exit
          if (msk(i,j) .eq. 0) r = 1
       end do
    end do
  end function amrex_mlndlap_any_zero


  subroutine amrex_mlndlap_set_dirichlet_mask (dmsk, dlo, dhi, omsk, olo, ohi, &
       domlo, domhi, bclo, bchi) bind(c,name='amrex_mlndlap_set_dirichlet_mask')
    integer, dimension(2) :: dlo, dhi, olo, ohi, domlo, domhi, bclo, bchi
    integer, intent(inout) :: dmsk(dlo(1):dhi(1),dlo(2):dhi(2))
    integer, intent(in   ) :: omsk(olo(1):ohi(1),olo(2):ohi(2))

    integer :: i,j
    
    do j = dlo(2), dhi(2)
       do i = dlo(1), dhi(1)
          if (any(omsk(i-1:i,j-1:j).eq.1)) then
             dmsk(i,j) = dirichlet
          else
             dmsk(i,j) = 0
          end if
       end do
    end do

    if (dlo(1) .eq. domlo(1)) then
       if (bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) then
          dmsk(dlo(1),:) = 0
       end if
    end if

    if (dhi(1) .eq. domhi(1)) then
       if (bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) then
          dmsk(dhi(1),:) = 0
       end if
    end if

    if (dlo(2) .eq. domlo(2)) then
       if (bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) then
          dmsk(:,dlo(2)) = 0
       end if
    end if

    if (dhi(2) .eq. domhi(2)) then
       if (bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) then
          dmsk(:,dhi(2)) = 0
       end if
    end if
    
    if (dlo(1) .eq. domlo(1)) then
       if (bclo(1) .eq. amrex_lo_dirichlet) then
          dmsk(dlo(1),:) = dirichlet
       end if
    end if

    if (dhi(1) .eq. domhi(1)) then
       if (bchi(1) .eq. amrex_lo_dirichlet) then
          dmsk(dhi(1),:) = dirichlet
       end if
    end if

    if (dlo(2) .eq. domlo(2)) then
       if (bclo(2) .eq. amrex_lo_dirichlet) then
          dmsk(:,dlo(2)) = dirichlet
       end if
    end if

    if (dhi(2) .eq. domhi(2)) then
       if (bchi(2) .eq. amrex_lo_dirichlet) then
          dmsk(:,dhi(2)) = dirichlet
       end if
    end if
    
  end subroutine amrex_mlndlap_set_dirichlet_mask


  subroutine amrex_mlndlap_divu_fine_contrib (clo, chi, lo, hi, rhs, rlo, rhi, &
       vel, vlo, vhi, frh, flo, fhi, msk, mlo, mhi, dxinv, ndlo, ndhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_divu_fine_contrib')
    integer, dimension(2), intent(in) :: clo, chi, lo, hi, rlo, rhi, vlo, vhi, flo, fhi, mlo, mhi, &
         ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1),vlo(2):vhi(2),2)
    real(amrex_real), intent(in   ) :: frh(flo(1):fhi(1),flo(2):fhi(2))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))

    integer :: i, j, ii, jj, side, id, jd, ic, jc
    real(amrex_real) :: facx, facy, facx_s, facy_s
    real(amrex_real), parameter :: rfd = 0.25d0
    real(amrex_real), parameter :: chip = 0.5d0
    real(amrex_real), parameter :: chip2 = 0.25d0

    ! note that dxinv is fine dxinv
    facx = 0.5d0*dxinv(1) * rfd
    facy = 0.5d0*dxinv(2) * rfd

    ! note that lo = 2*clo and hi = 2*chi

    do side = 0, 1
       if (side .eq. 0) then
          j = clo(2)
          jj = lo(2)
          facy_s = facy
          jd = jj+1
          jc = jj
       else
          j = chi(2)
          jj = hi(2)
          facy_s = -facy
          jd = jj-1
          jc = jj-1
       end if

       i = clo(1)
       ii = lo(1)
       if (msk(ii,jj) .eq. dirichlet) then
          rhs(i,j) = rhs(i,j) + facx*vel(ii,jc,1) + facy_s*vel(ii,jc,2) &
               + chip * (facx  *(vel(ii+1,jc,1)-vel(ii,jc,1)) &
               &       + facy_s*(vel(ii+1,jc,2)+vel(ii,jc,2))) &
               + rfd*chip2*frh(ii+1,jd)
       end if

       do i = clo(1)+1, chi(1)-1
          ii = 2*i
          if (msk(ii,jj) .eq. dirichlet) then
             rhs(i,j) = rhs(i,j) + facx  *(vel(ii,jc,1)-vel(ii-1,jc,1)) &
                  &              + facy_s*(vel(ii,jc,2)+vel(ii-1,jc,2)) &
                  + chip * (facx  *(vel(ii-1,jc,1)-vel(ii-2,jc,1)) &
                  &       + facy_s*(vel(ii-1,jc,2)+vel(ii-2,jc,2))) &
                  + chip * (facx  *(vel(ii+1,jc,1)-vel(ii,jc,1)) &
                  &       + facy_s*(vel(ii+1,jc,2)+vel(ii,jc,2))) &
                  + rfd*(chip2*frh(ii-1,jd) + chip*frh(ii,jd) + chip2*frh(ii+1,jd))
          end if
       end do

       i = chi(1)
       ii = hi(1)
       if (msk(ii,jj) .eq. dirichlet) then
          rhs(i,j) = rhs(i,j) - facx*vel(ii-1,jc,1) + facy_s*vel(ii-1,jc,2) &
               + chip * (facx  *(vel(ii-1,jc,1)-vel(ii-2,jc,1)) &
               &       + facy_s*(vel(ii-1,jc,2)+vel(ii-2,jc,2))) &
               + rfd*chip2*frh(ii-1,jd)
       end if
    end do

    do side = 0, 1
       if (side .eq. 0) then
          i = clo(1)
          ii = lo(1)
          facx_s = facx
          id = ii+1
          ic = ii
       else
          i = chi(1)
          ii = hi(1)
          facx_s = -facx
          id = ii-1
          ic = ii-1
       end if

       j = clo(2)
       jj = lo(2)
       if (msk(ii,jj) .eq. dirichlet) then
          rhs(i,j) = rhs(i,j) + chip * (facx_s*(vel(ic,jj+1,1)+vel(ic,jj,1)) &
               &                      + facy  *(vel(ic,jj+1,2)-vel(ic,jj,2)))
       end if

       do j = clo(2)+1, chi(2)-1
          jj = 2*j
          if (msk(ii,jj) .eq. dirichlet) then
             rhs(i,j) = rhs(i,j) + facx_s*(vel(ic,jj,1)+vel(ic,jj-1,1)) &
                  &              + facy  *(vel(ic,jj,2)-vel(ic,jj-1,2)) &
                  + chip * (facx_s*(vel(ic,jj-1,1)+vel(ic,jj-2,1)) &
                  &       + facy  *(vel(ic,jj-1,2)-vel(ic,jj-2,2))) &
                  + chip * (facx_s*(vel(ic,jj+1,1)+vel(ic,jj,1)) &
                  &       + facy  *(vel(ic,jj+1,2)-vel(ic,jj,2))) &
                  + rfd*(chip2*frh(id,jj-1) + chip*frh(id,jj) + chip2*frh(id,jj+1))
          end if
       end do

       j = chi(2)
       jj = hi(2)
       if (msk(ii,jj) .eq. dirichlet) then
          rhs(i,j) = rhs(i,j) + chip * (facx_s*(vel(ic,jj-1,1)+vel(ic,jj-2,1)) &
               &                      + facy  *(vel(ic,jj-1,2)-vel(ic,jj-2,2)))
       end if
    end do

    ! xxxxx what do we do at physical boundaries?

  end subroutine amrex_mlndlap_divu_fine_contrib


  subroutine amrex_mlndlap_divu_cf_contrib (lo, hi,  rhs, rlo, rhi, vel, vlo, vhi, dmsk, mlo, mhi, &
       fmsk, flo, fhi, fc, clo, chi, dxinv, ndlo, ndhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_divu_cf_contrib')
    integer, dimension(2), intent(in) :: lo, hi, rlo, rhi, vlo, vhi, mlo, mhi, flo, fhi, &
         clo, chi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1),vlo(2):vhi(2),2)
    real(amrex_real), intent(in   ) :: fc (clo(1):chi(1),clo(2):chi(2))
    integer, intent(in) :: dmsk(mlo(1):mhi(1),mlo(2):mhi(2))
    integer, intent(in) :: fmsk(flo(1):fhi(1),flo(2):fhi(2))

    integer :: i,j
    real(amrex_real) :: facx, facy

    facx = 0.5d0*dxinv(1)
    facy = 0.5d0*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (dmsk(i,j) .ne. dirichlet) then
             if (any(fmsk(i-1:i,j-1:j).eq.1) .and. any(fmsk(i-1:i,j-1:j).eq.0)) then
                rhs(i,j) = fc(i,j) &
                     + (1.d0-fmsk(i-1,j-1)) * (-facx*vel(i-1,j-1,1) - facy*vel(i-1,j-1,2)) &
                     + (1.d0-fmsk(i  ,j-1)) * ( facx*vel(i  ,j-1,1) - facy*vel(i  ,j-1,2)) &
                     + (1.d0-fmsk(i-1,j  )) * (-facx*vel(i-1,j  ,1) + facy*vel(i-1,j  ,2)) &
                     + (1.d0-fmsk(i  ,j  )) * ( facx*vel(i  ,j  ,1) + facy*vel(i  ,j  ,2))
             end if
          end if
       end do
    end do

    ! xxxxx what do we do at physical boundaries?
  end subroutine amrex_mlndlap_divu_cf_contrib


  subroutine amrex_mlndlap_res_fine_contrib (clo, chi, lo, hi, f, flo, fhi, x, xlo, xhi, &
       sig, slo, shi, res, rlo, rhi, rhs, hlo, hhi, msk, mlo, mhi, dxinv) &
       bind(c,name='amrex_mlndlap_res_fine_contrib')
    integer, dimension(2), intent(in) :: clo, chi, lo, hi, flo, fhi, xlo, xhi, slo, shi, &
         rlo, rhi, hlo, hhi, mlo, mhi
    real(amrex_real), intent(inout) :: f  (flo(1):fhi(1),flo(2):fhi(2))
    real(amrex_real), intent(in   ) :: x  (xlo(1):xhi(1),xlo(2):xhi(2))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1),slo(2):shi(2))
    real(amrex_real), intent(in   ) :: res(rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(in   ) :: rhs(hlo(1):hhi(1),hlo(2):hhi(2))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2))
    real(amrex_real), intent(in) :: dxinv(2)

    integer :: i, j, ii, jj, side, jd, jc, id, ic
    real(amrex_real) :: facx, facy
    real(amrex_real), parameter :: rfd = 0.25d0
    real(amrex_real), parameter :: chip = 0.5d0
    real(amrex_real), parameter :: chip2 = 0.25d0

    facx = (1.d0/6.d0)*dxinv(1)*dxinv(1) * rfd
    facy = (1.d0/6.d0)*dxinv(2)*dxinv(2) * rfd

    ! note that lo = 2*clo and hi = 2*chi

    do side = 0, 1
       if (side .eq. 0) then
          j = clo(2)
          jj = lo(2)
          jd = jj+1
          jc = jj
       else
          j = chi(2)
          jj = hi(2)
          jd = jj-1
          jc = jj-1
       end if

       i = clo(1)
       ii = lo(1)
       if (msk(ii,jj) .eq. dirichlet) then
          f(i,j) = f(i,j) + sig(ii,jc)*(facx*(2.d0*(x(ii+1,jj)-x(ii  ,jj))  &
               &                             +     (x(ii+1,jd)-x(ii  ,jd))) &
               &                      + facy*(2.d0*(x(ii  ,jd)-x(ii  ,jj))  &
               &                             +     (x(ii+1,jd)-x(ii+1,jj)))) &
               + chip * sig(ii  ,jc)*(facx*(2.d0*(x(ii  ,jj)-x(ii+1,jj))  &
               &                           +     (x(ii  ,jd)-x(ii+1,jd))) &
               &                    + facy*(2.d0*(x(ii+1,jd)-x(ii+1,jj))  &
               &                           +     (x(ii  ,jd)-x(ii  ,jj)))) &
               + chip * sig(ii+1,jc)*(facx*(2.d0*(x(ii+2,jj)-x(ii+1,jj))  &
               &                           +     (x(ii+2,jd)-x(ii+1,jd))) &
               &                    + facy*(2.d0*(x(ii+1,jd)-x(ii+1,jj))  &
               &                           +     (x(ii+2,jd)-x(ii+2,jj)))) &
               + rfd*chip2*(rhs(ii+1,jd)-res(ii+1,jd))
       end if

       do i = clo(1)+1, chi(1)-1
          ii = 2*i
          if (msk(ii,jj) .eq. dirichlet) then
             f(i,j) = f(i,j) + sig(ii-1,jc)*(facx*(2.d0*(x(ii-1,jj)-x(ii  ,jj))  &
                  &                               +     (x(ii-1,jd)-x(ii  ,jd))) &
                  &                        + facy*(2.d0*(x(ii  ,jd)-x(ii  ,jj))  &
                  &                               +     (x(ii-1,jd)-x(ii-1,jj)))) &
                  &          + sig(ii  ,jc)*(facx*(2.d0*(x(ii+1,jj)-x(ii  ,jj))  &
                  &                               +     (x(ii+1,jd)-x(ii  ,jd))) &
                  &                        + facy*(2.d0*(x(ii  ,jd)-x(ii  ,jj))  &
                  &                               +     (x(ii+1,jd)-x(ii+1,jj)))) &
                  + chip * sig(ii-2,jc)*(facx*(2.d0*(x(ii-2,jj)-x(ii-1,jj))  &
                  &                           +     (x(ii-2,jd)-x(ii-1,jd))) &
                  &                    + facy*(2.d0*(x(ii-1,jd)-x(ii-1,jj))  &
                  &                           +     (x(ii-2,jd)-x(ii-2,jj)))) &
                  + chip * sig(ii-1,jc)*(facx*(2.d0*(x(ii  ,jj)-x(ii-1,jj))  &
                  &                           +     (x(ii  ,jd)-x(ii-1,jd))) &
                  &                    + facy*(2.d0*(x(ii-1,jd)-x(ii-1,jj))  &
                  &                           +     (x(ii  ,jd)-x(ii  ,jj)))) &
                  + chip * sig(ii  ,jc)*(facx*(2.d0*(x(ii  ,jj)-x(ii+1,jj))  &
                  &                           +     (x(ii  ,jd)-x(ii+1,jd))) &
                  &                    + facy*(2.d0*(x(ii+1,jd)-x(ii+1,jj))  &
                  &                           +     (x(ii  ,jd)-x(ii  ,jj)))) &
                  + chip * sig(ii+1,jc)*(facx*(2.d0*(x(ii+2,jj)-x(ii+1,jj))  &
                  &                           +     (x(ii+2,jd)-x(ii+1,jd))) &
                  &                    + facy*(2.d0*(x(ii+1,jd)-x(ii+1,jj))  &
                  &                           +     (x(ii+2,jd)-x(ii+2,jj)))) &
                  + rfd*(chip2*(rhs(ii-1,jd)-res(ii-1,jd))  &
                  &     +chip *(rhs(ii  ,jd)-res(ii  ,jd))  &
                  &     +chip2*(rhs(ii+1,jd)-res(ii+1,jd)))
          end if
       end do

       i = chi(1)
       ii = hi(1)
       if (msk(ii,jj) .eq. dirichlet) then
          f(i,j) = f(i,j) + sig(ii-1,jc)*(facx*(2.d0*(x(ii-1,jj)-x(ii  ,jj))  &
               &                               +     (x(ii-1,jd)-x(ii  ,jd))) &
               &                        + facy*(2.d0*(x(ii  ,jd)-x(ii  ,jj))  &
               &                               +     (x(ii-1,jd)-x(ii-1,jj)))) &
               + chip * sig(ii-2,jc)*(facx*(2.d0*(x(ii-2,jj)-x(ii-1,jj))  &
               &                           +     (x(ii-2,jd)-x(ii-1,jd))) &
               &                    + facy*(2.d0*(x(ii-1,jd)-x(ii-1,jj))  &
               &                           +     (x(ii-2,jd)-x(ii-2,jj)))) &
               + chip * sig(ii-1,jc)*(facx*(2.d0*(x(ii  ,jj)-x(ii-1,jj))  &
               &                           +     (x(ii  ,jd)-x(ii-1,jd))) &
               &                    + facy*(2.d0*(x(ii-1,jd)-x(ii-1,jj))  &
               &                           +     (x(ii  ,jd)-x(ii  ,jj)))) &
               + rfd*chip2*(rhs(ii-1,jd)-res(ii-1,jd))
       end if
    end do

    do side = 0, 1
       if (side .eq. 0) then
          i = clo(1)
          ii = lo(1)
          id = ii+1
          ic = ii
       else
          i = chi(1)
          ii = hi(1)
          id = ii-1
          ic = ii-1
       end if

       j = clo(2)
       jj = lo(2)
       if (msk(ii,jj) .eq. dirichlet) then
          f(i,j) = f(i,j)  &
               + chip * sig(ic,jj  )*(facx*(2.d0*(x(id,jj+1)-x(ii,jj+1))  &
               &                           +     (x(id,jj  )-x(ii,jj  ))) &
               &                    + facy*(2.d0*(x(ii,jj  )-x(ii,jj+1))  &
               &                           +     (x(id,jj  )-x(id,jj+1)))) &
               + chip * sig(ic,jj+1)*(facx*(2.d0*(x(id,jj+1)-x(ii,jj+1))  &
               &                           +     (x(id,jj+2)-x(ii,jj+2))) &
               &                    + facy*(2.d0*(x(ii,jj+2)-x(ii,jj+1))  &
               &                           +     (x(id,jj+2)-x(id,jj+1))))
       end if

       do j = clo(2)+1, chi(2)-1
          jj = 2*j
          if (msk(ii,jj) .eq. dirichlet) then
             f(i,j) = f(i,j) + sig(ic,jj-1)*(facx*(2.d0*(x(id,jj  )-x(ii,jj  ))  &
                  &                               +     (x(id,jj-1)-x(ii,jj-1))) &
                  &                        + facy*(2.d0*(x(ii,jj-1)-x(ii,jj  ))  &
                  &                               +     (x(id,jj-1)-x(id,jj  )))) &
                  &          + sig(ic,jj  )*(facx*(2.d0*(x(id,jj  )-x(ii,jj  ))  &
                  &                               +     (x(id,jj+1)-x(ii,jj+1))) &
                  &                        + facy*(2.d0*(x(ii,jj+1)-x(ii,jj  ))  &
                  &                               +     (x(id,jj+1)-x(id,jj  )))) &
                  + chip * sig(ic,jj-2)*(facx*(2.d0*(x(id,jj-1)-x(ii,jj-1))  &
                  &                           +     (x(id,jj-2)-x(ii,jj-2))) &
                  &                    + facy*(2.d0*(x(ii,jj-2)-x(ii,jj-1))  &
                  &                           +     (x(id,jj-2)-x(id,jj-1)))) &
                  + chip * sig(ic,jj-1)*(facx*(2.d0*(x(id,jj-1)-x(ii,jj-1))  &
                  &                           +     (x(id,jj  )-x(ii,jj  ))) &
                  &                    + facy*(2.d0*(x(ii,jj  )-x(ii,jj-1))  &
                  &                           +     (x(id,jj  )-x(id,jj-1)))) &
                  + chip * sig(ic,jj  )*(facx*(2.d0*(x(id,jj+1)-x(ii,jj+1))  &
                  &                           +     (x(id,jj  )-x(ii,jj  ))) &
                  &                    + facy*(2.d0*(x(ii,jj  )-x(ii,jj+1))  &
                  &                           +     (x(id,jj  )-x(id,jj+1)))) &
                  + chip * sig(ic,jj+1)*(facx*(2.d0*(x(id,jj+1)-x(ii,jj+1))  &
                  &                           +     (x(id,jj+2)-x(ii,jj+2))) &
                  &                    + facy*(2.d0*(x(ii,jj+2)-x(ii,jj+1))  &
                  &                           +     (x(id,jj+2)-x(id,jj+1)))) &
                  + rfd*(chip2*(rhs(id,jj-1)-res(id,jj-1)) &
                  &     +chip *(rhs(id,jj  )-res(id,jj  )) &
                  &     +chip2*(rhs(id,jj+1)-res(id,jj+1)))
          end if
       end do

       j = chi(2)
       jj = hi(2)
       if (msk(ii,jj) .eq. dirichlet) then
          f(i,j) = f(i,j) &
               + chip*sig(ic,jj-2)*(facx*(2.d0*(x(id,jj-1)-x(ii,jj-1))  &
               &                         +     (x(id,jj-2)-x(ii,jj-2))) &
               &                  + facy*(2.d0*(x(ii,jj-2)-x(ii,jj-1))  &
               &                         +     (x(id,jj-2)-x(id,jj-1)))) &
               + chip*sig(ic,jj-1)*(facx*(2.d0*(x(id,jj-1)-x(ii,jj-1))  &
               &                         +     (x(id,jj  )-x(ii,jj  ))) &
               &                  + facy*(2.d0*(x(ii,jj  )-x(ii,jj-1))  &
               &                         +     (x(id,jj  )-x(id,jj-1))))
       end if
    end do

  end subroutine amrex_mlndlap_res_fine_contrib


  subroutine amrex_mlndlap_res_cf_contrib (lo, hi, res, rlo, rhi, phi, phlo, phhi, &
       rhs, rhlo, rhhi, sig, slo, shi, dmsk, mlo, mhi, fmsk, flo, fhi, &
       fc, clo, chi, dxinv) &
       bind(c,name='amrex_mlndlap_res_cf_contrib')
    integer, dimension(2), intent(in) :: lo, hi, rlo, rhi, phlo, phhi, rhlo, rhhi, slo, shi, &
         mlo, mhi, flo, fhi, clo, chi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: res( rlo(1): rhi(1), rlo(2): rhi(2))
    real(amrex_real), intent(in   ) :: phi(phlo(1):phhi(1),phlo(2):phhi(2))
    real(amrex_real), intent(in   ) :: rhs(rhlo(1):rhhi(1),rhlo(2):rhhi(2))
    real(amrex_real), intent(in   ) :: sig( slo(1): shi(1), slo(2): shi(2))
    real(amrex_real), intent(inout) :: fc ( clo(1): chi(1), clo(2): chi(2))
    integer, intent(in) :: dmsk(mlo(1):mhi(1),mlo(2):mhi(2))
    integer, intent(in) :: fmsk(flo(1):fhi(1),flo(2):fhi(2))

    integer :: i,j
    real(amrex_real) :: Ax, facx, facy

    facx = (1.d0/6.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/6.d0)*dxinv(2)*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (dmsk(i,j) .ne. dirichlet) then
             if (any(fmsk(i-1:i,j-1:j).eq.0) .and. any(fmsk(i-1:i,j-1:j).eq.1)) then
                Ax = 0.d0
                if (fmsk(i-1,j-1) .eq. 0) then
                   Ax = Ax + sig(i-1,j-1)*(facx*(2.d0*(phi(i-1,j  )-phi(i  ,j  )) &
                        &                       +     (phi(i-1,j-1)-phi(i  ,j-1))) &
                        &                + facy*(2.d0*(phi(i  ,j-1)-phi(i  ,j  )) &
                        &                       +     (phi(i-1,j-1)-phi(i-1,j  ))))
                end if
                if (fmsk(i,j-1) .eq. 0) then
                   Ax = Ax + sig(i,j-1)*(facx*(2.d0*(phi(i+1,j  )-phi(i  ,j  )) &
                        &                     +     (phi(i+1,j-1)-phi(i  ,j-1))) &
                        &              + facy*(2.d0*(phi(i  ,j-1)-phi(i  ,j  )) &
                        &                     +     (phi(i+1,j-1)-phi(i+1,j  ))))
                end if
                if (fmsk(i-1,j) .eq. 0) then
                   Ax = Ax + sig(i-1,j)*(facx*(2.d0*(phi(i-1,j  )-phi(i  ,j  )) &
                        &                     +     (phi(i-1,j+1)-phi(i  ,j+1))) &
                        &              + facy*(2.d0*(phi(i  ,j+1)-phi(i  ,j  )) &
                        &                     +     (phi(i-1,j+1)-phi(i-1,j  ))))
                end if
                if (fmsk(i,j) .eq. 0) then
                   Ax = Ax + sig(i,j)*(facx*(2.d0*(phi(i+1,j  )-phi(i  ,j  )) &
                        &                  +      (phi(i+1,j+1)-phi(i  ,j+1))) &
                        &            + facy*(2.d0*(phi(i  ,j+1)-phi(i  ,j  )) &
                        &                  +      (phi(i+1,j+1)-phi(i+1,j  ))))
                end if
                res(i,j) = rhs(i,j) - Ax - fc(i,j)
             end if
          end if
       end do
    end do
  end subroutine amrex_mlndlap_res_cf_contrib


  subroutine amrex_mlndlap_fixup_res_mask (lo, hi, rmsk, rlo, rhi, fmsk, flo, fhi) &
       bind(c,name='amrex_mlndlap_fixup_res_mask')
    integer, dimension(2), intent(in) :: lo, hi, rlo, rhi, flo, fhi
    integer, intent(inout) :: rmsk(rlo(1):rhi(1),rlo(2):rhi(2))
    integer, intent(in   ) :: fmsk(flo(1):fhi(1),flo(2):fhi(2))

    integer :: i,j
    
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (any(fmsk(i-1:i,j-1:j).eq.0) .and. any(fmsk(i-1:i,j-1:j).eq.1)) then
             rmsk(i,j) = 1
          end if
       end do
    end do
  end subroutine amrex_mlndlap_fixup_res_mask

end module amrex_mlnodelap_2d_module
