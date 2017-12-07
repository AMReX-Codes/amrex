module amrex_mlnodelap_2d_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  use amrex_lo_bctypes_module, only : amrex_lo_dirichlet, amrex_lo_neumann, amrex_lo_inflow
  implicit none

  private
  public :: amrex_mlndlap_avgdown_coeff, amrex_mlndlap_fillbc_coeff, amrex_mlndlap_divu, &
       amrex_mlndlap_applybc, amrex_mlndlap_adotx, amrex_mlndlap_jacobi, &
       amrex_mlndlap_restriction, amrex_mlndlap_interpolation

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


  subroutine amrex_mlndlap_fillbc_coeff (sigma, slo, shi, dlo, dhi) &
       bind(c, name='amrex_mlndlap_fillbc_coeff')
    integer, dimension(2), intent(in) :: slo, shi, dlo, dhi
    real(amrex_real), intent(inout) :: sigma(slo(1):shi(1),slo(2):shi(2))

    integer :: ilo, ihi, jlo, jhi

    ilo = max(dlo(1), slo(1))
    ihi = min(dhi(1), shi(1))
    jlo = max(dlo(2), slo(2))
    jhi = min(dhi(2), shi(2))

    if (slo(1) .lt. dlo(1)) then
       sigma(dlo(1)-1,jlo:jhi) = sigma(dlo(1),jlo:jhi)
    end if
    
    if (shi(1) .gt. dhi(1)) then
       sigma(dhi(1)+1,jlo:jhi) = sigma(dhi(1),jlo:jhi)
    end if

    if (slo(2) .lt. dlo(2)) then
       sigma(ilo:ihi,dlo(2)-1) = sigma(ilo:ihi,dlo(2))
    end if

    if (shi(2) .gt. dhi(2)) then
       sigma(ilo:ihi,dhi(2)+1) = sigma(ilo:ihi,dhi(2))
    end if

    if (slo(1) .lt. dlo(1) .and. slo(2) .lt. dlo(2)) then
       sigma(dlo(1)-1,dlo(2)-1) = sigma(dlo(1),dlo(2))
    end if

    if (shi(1) .gt. dhi(1) .and. slo(2) .lt. dlo(2)) then
       sigma(dhi(1)+1,dlo(2)-1) = sigma(dhi(1),dlo(2))
    end if

    if (slo(1) .lt. dlo(1) .and. shi(2) .gt. dhi(2)) then
       sigma(dlo(1)-1,dhi(2)+1) = sigma(dlo(1),dhi(2))
    end if

    if (shi(1) .gt. dhi(1) .and. shi(2) .gt. dhi(2)) then
       sigma(dhi(1)+1,dhi(2)+1) = sigma(dhi(1),dhi(2))
    end if

  end subroutine amrex_mlndlap_fillbc_coeff


  subroutine amrex_mlndlap_divu (lo, hi, rhs, rlo, rhi, vel, vlo, vhi, dxinv, ndlo, ndhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_divu')
    integer, dimension(2), intent(in) :: lo, hi, rlo, rhi, vlo, vhi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1),vlo(2):vhi(2),2)

    integer :: i,j
    real(amrex_real) :: facx, facy

    facx = 0.5d0*dxinv(1)
    facy = 0.5d0*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          rhs(i,j) = facx*(-vel(i-1,j-1,1)+vel(i,j-1,1)-vel(i-1,j,1)+vel(i,j,1)) &
               &   + facy*(-vel(i-1,j-1,2)-vel(i,j-1,2)+vel(i-1,j,2)+vel(i,j,2))
       end do
    end do

    if (lo(1) .eq. ndlo(1)) then
       if (bclo(1) .eq. amrex_lo_dirichlet) then
          rhs(lo(1),lo(2):hi(2)) = 0.d0
       else if (bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) then 
          rhs(lo(1),lo(2):hi(2)) = 2.d0*rhs(lo(1),lo(2):hi(2))
       end if
    end if

    if (hi(1) .eq. ndhi(1)) then
       if (bchi(1) .eq. amrex_lo_dirichlet) then
          rhs(hi(1),lo(2):hi(2)) = 0.d0
       else if (bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) then
          rhs(hi(1),lo(2):hi(2)) = 2.d0*rhs(hi(1),lo(2):hi(2))
       end if
    end if

    if (lo(2) .eq. ndlo(2)) then
       if (bclo(2) .eq. amrex_lo_dirichlet) then
          rhs(lo(1):hi(1),lo(2)) = 0.d0
       else if (bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) then
          rhs(lo(1):hi(1),lo(2)) = 2.d0*rhs(lo(1):hi(1),lo(2))
       end if
    end if

    if (hi(2) .eq. ndhi(2)) then
       if (bchi(2) .eq. amrex_lo_dirichlet) then
          rhs(lo(1):hi(1),hi(2)) = 0.d0
       else if (bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) then
          rhs(lo(1):hi(1),hi(2)) = 2.d0*rhs(lo(1):hi(1),hi(2))
       end if
    end if

  end subroutine amrex_mlndlap_divu


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

    ! dirichlet

    if (bclo(1) .eq. amrex_lo_dirichlet .and. hlo(1) .lt. dlo(1)) then
       phi(dlo(1)-1:dlo(1),jlo:jhi) = 0.d0
    end if

    if (bchi(1) .eq. amrex_lo_dirichlet .and. hhi(1) .gt. dhi(1)) then
       phi(dhi(1):dhi(1)+1,jlo:jhi) = 0.d0
    end if

    if (bclo(2) .eq. amrex_lo_dirichlet .and. hlo(2) .lt. dlo(2)) then
       phi(ilo:ihi,dlo(2)-1:dlo(2)) = 0.d0
    end if

    if (bchi(2) .eq. amrex_lo_dirichlet .and. hhi(2) .gt. dhi(2)) then
       phi(ilo:ihi,dhi(2)+1) = phi(ilo:ihi,dhi(2)-1) 
    end if

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
       if (bclo(1) .eq. amrex_lo_dirichlet .or. bclo(2) .eq. amrex_lo_dirichlet) then
          phi(dlo(1)-1,dlo(2)-1) = 0.d0
       else
          phi(dlo(1)-1,dlo(2)-1) = phi(dlo(1)+1,dlo(2)+1)
       end if
    end if

    if (hhi(1) .gt. dhi(1) .and. hlo(2) .lt. dlo(2)) then
       if (bchi(1) .eq. amrex_lo_dirichlet .or. bclo(2) .eq. amrex_lo_dirichlet) then
          phi(dhi(1)+1,dlo(2)-1) = 0.d0
       else
          phi(dhi(1)+1,dlo(2)-1) = phi(dhi(1)-1,dlo(2)+1)
       end if
    end if

    if (hlo(1) .lt. dlo(1) .and. hhi(2) .gt. dhi(2)) then
       if (bclo(1) .eq. amrex_lo_dirichlet .or. bchi(2) .eq. amrex_lo_dirichlet) then
          phi(dlo(1)-1,dhi(2)+1) = 0.d0
       else
          phi(dlo(1)-1,dhi(2)+1) = phi(dlo(1)+1,dhi(2)-1)
       end if
    end if

    if (hhi(1) .gt. dhi(1) .and. hhi(2) .gt. dhi(2)) then
       if (bchi(1) .eq. amrex_lo_dirichlet .or. bchi(2) .eq. amrex_lo_dirichlet) then
          phi(dhi(1)+1,dhi(2)+1) = 0.d0
       else
          phi(dhi(1)+1,dhi(2)+1) = phi(dhi(1)-1,dhi(2)-1)
       end if
    end if

  end subroutine amrex_mlndlap_applybc


  subroutine amrex_mlndlap_adotx (lo, hi, y, ylo, yhi, x, xlo, xhi, &
       sx, sxlo, sxhi, sy, sylo, syhi, dg, dlo, dhi, dxinv) bind(c,name='amrex_mlndlap_adotx')
    integer, dimension(2), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, sxlo, sxhi, sylo, syhi, dlo, dhi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) ::  y( ylo(1): yhi(1), ylo(2): yhi(2))
    real(amrex_real), intent(in   ) ::  x( xlo(1): xhi(1), xlo(2): xhi(2))
    real(amrex_real), intent(in   ) :: sx(sxlo(1):sxhi(1),sxlo(2):sxhi(2))
    real(amrex_real), intent(in   ) :: sy(sylo(1):syhi(1),sylo(2):syhi(2))
    real(amrex_real)                :: dg( dlo(1): dhi(1), dlo(2): dhi(2), 2)
    
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
          y(i,j) = facx*(-sx(i-1,j-1)*(dg(i-1,j-1,1)+2.d0*dg(i-1,j  ,1)) &
               &         +sx(i  ,j-1)*(dg(i  ,j-1,1)+2.d0*dg(i  ,j  ,1)) &
               &         -sx(i-1,j  )*(dg(i-1,j+1,1)+2.d0*dg(i-1,j  ,1)) &
               &         +sx(i  ,j  )*(dg(i  ,j+1,1)+2.d0*dg(i  ,j  ,1))) &
               +   facy*(-sy(i-1,j-1)*(dg(i-1,j-1,2)+2.d0*dg(i  ,j-1,2)) &
               &         -sy(i  ,j-1)*(dg(i+1,j-1,2)+2.d0*dg(i  ,j-1,2)) &
               &         +sy(i-1,j  )*(dg(i-1,j  ,2)+2.d0*dg(i  ,j  ,2)) &
               &         +sy(i  ,j  )*(dg(i+1,j  ,2)+2.d0*dg(i  ,j  ,2)))
       end do
    end do

  end subroutine amrex_mlndlap_adotx  


  subroutine amrex_mlndlap_jacobi (lo, hi, sol, slo, shi, Ax, alo, ahi, rhs, rlo, rhi, &
       sx, sxlo, sxhi, sy, sylo, syhi, dxinv) bind(c,name='amrex_mlndlap_jacobi')
    integer, dimension(2),intent(in) :: lo,hi,slo,shi,alo,ahi,rlo,rhi,sxlo,sxhi,sylo,syhi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: sol( slo(1): shi(1), slo(2): shi(2))
    real(amrex_real), intent(in   ) :: Ax ( alo(1): ahi(1), alo(2): ahi(2))
    real(amrex_real), intent(in   ) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2))
    real(amrex_real), intent(in   ) :: sx (sxlo(1):sxhi(1),sxlo(2):sxhi(2))
    real(amrex_real), intent(in   ) :: sy (sylo(1):syhi(1),sylo(2):syhi(2))

    integer :: i,j
    real(amrex_real) :: facx, facy
    real(amrex_real), parameter :: omega = 2.d0/3.d0

    facx = -2.d0 * (1.d0/6.d0)*dxinv(1)*dxinv(1)
    facy = -2.d0 * (1.d0/6.d0)*dxinv(2)*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          sol(i,j) = (1.d0-omega)*sol(i,j) + omega * (rhs(i,j) - Ax(i,j)) &
               / (facx*(sx(i-1,j-1)+sx(i,j-1)+sx(i-1,j)+sx(i,j)) &
               +  facy*(sy(i-1,j-1)+sy(i,j-1)+sy(i-1,j)+sy(i,j)))
       end do
    end do
  end subroutine amrex_mlndlap_jacobi


  subroutine amrex_mlndlap_restriction (lo, hi, crse, clo, chi, fine, flo, fhi, dlo, dhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_restriction')
    integer, dimension(2), intent(in) :: lo, hi, clo, chi, flo, fhi, dlo, dhi, bclo, bchi
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2))
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1),flo(2):fhi(2))

    integer :: i,j, ii, jj
    real(amrex_real), parameter :: fac = 1.d0/16.d0

    do    j = lo(2), hi(2)
       jj = 2*j
       do i = lo(1), hi(1)
          ii = 2*i
          crse(i,j) = fac*(fine(ii-1,jj-1) + 2.d0*fine(ii  ,jj-1) +      fine(ii+1,jj-1) &
               +      2.d0*fine(ii-1,jj  ) + 4.d0*fine(ii  ,jj  ) + 2.d0*fine(ii+1,jj  ) &
               +           fine(ii-1,jj+1) + 2.d0*fine(ii  ,jj+1) +      fine(ii+1,jj+1))
       end do
    end do

    if (bclo(1) .eq. amrex_lo_dirichlet .and. lo(1) .eq. dlo(1)) then
       crse(lo(1),lo(2):hi(2)) = 0.d0
    end if

    if (bchi(1) .eq. amrex_lo_dirichlet .and. hi(1) .eq. dhi(1)) then
       crse(hi(1),lo(2):hi(2)) = 0.d0
    end if

    if (bclo(2) .eq. amrex_lo_dirichlet .and. lo(2) .eq. dlo(2)) then
       crse(lo(1):hi(1),lo(2)) = 0.d0
    end if

    if (bchi(2) .eq. amrex_lo_dirichlet .and. hi(2) .eq. dhi(2)) then
       crse(lo(1):hi(1),hi(2)) = 0.d0
    end if
  end subroutine amrex_mlndlap_restriction


  subroutine amrex_mlndlap_interpolation (clo, chi, fine, fflo, ffhi, crse, cflo, cfhi, &
       sigx, sxlo, sxhi, sigy, sylo, syhi) bind(c,name='amrex_mlndlap_interpolation')
    integer, dimension(2), intent(in) :: clo,chi,fflo,ffhi,cflo,cfhi,sxlo,sxhi,sylo,syhi
    real(amrex_real), intent(in   ) :: crse(cflo(1):cfhi(1),cflo(2):cfhi(2))
    real(amrex_real), intent(inout) :: fine(fflo(1):ffhi(1),fflo(2):ffhi(2))
    real(amrex_real), intent(in   ) :: sigx(sxlo(1):sxhi(1),sxlo(2):sxhi(2))
    real(amrex_real), intent(in   ) :: sigy(sylo(1):syhi(1),sylo(2):syhi(2))

    integer :: flo(2), fhi(2), i,j
    logical :: interpx
    real(amrex_real) :: wxm, wxp, wym, wyp

    flo = 2*clo
    fhi = 2*chi

    do    j = clo(2), chi(2)
       do i = clo(1), chi(1)
          fine(2*i,2*j) = crse(i,j)
       end do
    end do

    interpx = .false.
    do j = flo(2), fhi(2)
       interpx = .not.interpx
       if (interpx) then
          ! interp in x-direction
          do i = flo(1)+1, fhi(1), 2
             wxm = sigx(i-1,j-1) + sigx(i-1,j)
             wxp = sigx(i  ,j-1) + sigx(i  ,j)
             fine(i,j) = (wxm*fine(i-1,j) + wxp*fine(i+1,j)) / (wxm+wxp)
          end do
       else
          ! interp in y-direction
          do i = flo(1), fhi(1), 2
             wym = sigy(i-1,j-1) + sigy(i,j-1)
             wyp = sigy(i-1,j  ) + sigy(i,j  )
             fine(i,j) = (wym*fine(i,j-1) + wyp*fine(i,j+1)) / (wym+wyp)
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

  end subroutine amrex_mlndlap_interpolation

end module amrex_mlnodelap_2d_module
