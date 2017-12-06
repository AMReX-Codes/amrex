module amrex_mlnodelap_2d_module

  use amrex_fort_module, only : amrex_real
  use amrex_lo_bctypes_module, only : amrex_lo_dirichlet, amrex_lo_neumann, amrex_lo_inflow
  implicit none

  private
  public :: amrex_mlndlap_avgdown_coeff, amrex_mlndlap_divu, &
       amrex_mlndlap_adotx

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

end module amrex_mlnodelap_2d_module
