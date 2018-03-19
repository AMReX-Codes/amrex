
#include "AMReX_LO_BCTYPES.H"

module amrex_mllinop_nd_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real, amrex_spacedim
  implicit none

#if (AMREX_SPACEDIM == 1)
  integer, parameter :: xlo_dir = 0
  integer, parameter :: xhi_dir = 1
#elif (AMREX_SPACEDIM == 2)
  integer, parameter :: xlo_dir = 0
  integer, parameter :: ylo_dir = 1
  integer, parameter :: xhi_dir = 2
  integer, parameter :: yhi_dir = 3
#else
  integer, parameter :: xlo_dir = 0
  integer, parameter :: ylo_dir = 1
  integer, parameter :: zlo_dir = 2
  integer, parameter :: xhi_dir = 3
  integer, parameter :: yhi_dir = 4
  integer, parameter :: zhi_dir = 5
#endif
  
  private
  public :: amrex_mllinop_apply_bc, amrex_mllinop_comp_interp_coef0, amrex_mllinop_apply_metric

contains

  subroutine amrex_mllinop_apply_bc (lo, hi, phi, hlo, hhi, mask, mlo, mhi, &
       cdir, bct, bcl, bcval, blo, bhi, maxorder, dxinv, inhomog, nc, cross) &
       bind(c,name='amrex_mllinop_apply_bc')
    integer, dimension(3), intent(in) :: lo, hi, hlo, hhi, mlo, mhi, blo, bhi
    integer, value, intent(in) :: cdir, bct, maxorder, inhomog, nc, cross
    real(amrex_real), value, intent(in) :: bcl
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) ::  phi (hlo(1):hhi(1),hlo(2):hhi(2),hlo(3):hhi(3),nc)
    integer         , intent(in   ) :: mask (mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    real(amrex_real), intent(in   ) :: bcval(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),nc)

    integer :: i, j, k, idim, lenx, m
    logical :: inhomogeneous
    real(amrex_real) ::    x(-1:maxorder-2)
    real(amrex_real) :: coef(-1:maxorder-2), coef2(-maxorder+2:1)
    real(amrex_real), parameter :: xInt = -0.5D0
    real(amrex_real) :: fac

    integer :: n
    

    inhomogeneous = (inhomog .ne. 0)

    do n = 1, nc
       if (bct == LO_NEUMANN .or. bct == LO_REFLECT_ODD) then

          if (bct == LO_NEUMANN) then
             fac = 1.d0
          else
             fac = -1.d0
          end if

          select case (cdir)
          case (xlo_dir)
             do    k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   if (mask(lo(1)-1,j,k) .gt. 0) then
                      phi(lo(1)-1,j,k,n) = fac*phi(lo(1),j,k,n)
                   end if
                end do
             end do
          case (xhi_dir)
             do    k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   if (mask(hi(1)+1,j,k) .gt. 0) then
                      phi(hi(1)+1,j,k,n) = fac*phi(hi(1),j,k,n)
                   end if
                end do
             end do
#if (AMREX_SPACEDIM >= 2)
          case (ylo_dir)  
             do    k = lo(3), hi(3)
                do i = lo(1), hi(1)
                   if (mask(i,lo(2)-1,k) .gt. 0) then
                      phi(i,lo(2)-1,k,n) = fac*phi(i,lo(2),k,n)
                   end if
                end do
             end do
          case (yhi_dir)
             do    k = lo(3), hi(3)
                do i = lo(1), hi(1)
                   if (mask(i,hi(2)+1,k) .gt. 0) then
                      phi(i,hi(2)+1,k,n) = fac*phi(i,hi(2),k,n)
                   end if
                end do
             end do
#if (AMREX_SPACEDIM == 3)
          case (zlo_dir)
             do    j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if (mask(i,j,lo(3)-1) .gt. 0) then
                      phi(i,j,lo(3)-1,n) = fac*phi(i,j,lo(3),n)
                   end if
                end do
             end do
          case (zhi_dir)
             do    j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if (mask(i,j,hi(3)+1) .gt. 0) then
                      phi(i,j,hi(3)+1,n) = fac*phi(i,j,hi(3),n)
                   end if
                end do
             end do
#endif
#endif
          end select

       else if (bct == LO_DIRICHLET) then

          idim = mod(cdir,amrex_spacedim) + 1 ! cdir starts with 0; idim starts with 1
          lenx = MIN(hi(idim)-lo(idim), maxorder-2)

          x(-1) = -bcl*dxinv(idim)
          do m=0,maxorder-2
             x(m) = m + 0.5D0
          end do

          call polyInterpCoeff(xInt, x, lenx+2, coef)
          do m = -lenx, 1
             coef2(m) = coef(-m)
          end do

          select case (cdir)
          case (xlo_dir)
             do    k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   if (mask(lo(1)-1,j,k) .gt. 0) then
                      phi(lo(1)-1,j,k,n) = sum(phi(lo(1):lo(1)+lenx,j,k,n)*coef(0:lenx))
                      if (inhomogeneous) then
                         phi(lo(1)-1,j,k,n) = phi(lo(1)-1,j,k,n) + bcval(lo(1)-1,j,k,n)*coef(-1)
                      end if
                   end if
                end do
             end do
          case (xhi_dir)
             do    k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   if (mask(hi(1)+1,j,k) .gt. 0) then
                      phi(hi(1)+1,j,k,n) = sum(phi(hi(1)-lenx:hi(1),j,k,n)*coef2(-lenx:0))
                      if (inhomogeneous) then
                         phi(hi(1)+1,j,k,n) = phi(hi(1)+1,j,k,n) + bcval(hi(1)+1,j,k,n)*coef2(1)
                      end if
                   end if
                end do
             end do
#if (AMREX_SPACEDIM >= 2)
          case (ylo_dir)
             do    k = lo(3), hi(3)
                do i = lo(1), hi(1)
                   if (mask(i,lo(2)-1,k) .gt. 0) then
                      phi(i,lo(2)-1,k,n) = sum(phi(i,lo(2):lo(2)+lenx,k,n)*coef(0:lenx))
                      if (inhomogeneous) then
                         phi(i,lo(2)-1,k,n) = phi(i,lo(2)-1,k,n) + bcval(i,lo(2)-1,k,n)*coef(-1)
                      end if
                   end if
                end do
             end do
          case (yhi_dir)
             do    k = lo(3), hi(3)
                do i = lo(1), hi(1)
                   if (mask(i,hi(2)+1,k) .gt. 0) then
                      phi(i,hi(2)+1,k,n) = sum(phi(i,hi(2)-lenx:hi(2),k,n)*coef2(-lenx:0))
                      if (inhomogeneous) then
                         phi(i,hi(2)+1,k,n) = phi(i,hi(2)+1,k,n) + bcval(i,hi(2)+1,k,n)*coef2(1)
                      end if
                   end if
                end do
             end do
#if (AMREX_SPACEDIM == 3)
          case (zlo_dir)
             do    j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if (mask(i,j,lo(3)-1) .gt. 0) then
                      phi(i,j,lo(3)-1,n) = sum(phi(i,j,lo(3):lo(3)+lenx,n)*coef(0:lenx))
                      if (inhomogeneous) then
                         phi(i,j,lo(3)-1,n) = phi(i,j,lo(3)-1,n) + bcval(i,j,lo(3)-1,n)*coef(-1)
                      end if
                   end if
                end do
             end do
          case (zhi_dir)
             do    j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if (mask(i,j,hi(3)+1) .gt. 0) then
                      phi(i,j,hi(3)+1,n) = sum(phi(i,j,hi(3)-lenx:hi(3),n)*coef2(-lenx:0))
                      if (inhomogeneous) then
                         phi(i,j,hi(3)+1,n) = phi(i,j,hi(3)+1,n) + bcval(i,j,hi(3)+1,n)*coef2(1)
                      end if
                   end if
                end do
             end do
#endif
#endif
          end select

       else
          call amrex_error("amrex_mllinop_nd_module: unknown bc");
       end if

       ! Fill corners with averages for non-cross stencil
#if (AMREX_SPACEDIM == 2)
       if (cross .eq. 0) then
          ! The iteration over faces is always in the order of xlo, xhi, ylo and yhi.
          ! No need to do anything for xlo and xhi, because at that time, ylo and yhi
          ! have not been filled.
          select case (cdir)
          case (ylo_dir)
             do k = lo(3), hi(3)
                if (mask(lo(1)-1,lo(2)-1,k) .gt. 0) then
                   phi(lo(1)-1,lo(2)-1,k,n) = & ! Southwest
                        0.5d0*(2.d0*phi(lo(1),lo(2)-1,k,n) - phi(lo(1)+1,lo(2)-1,k,n)) + &
                        0.5d0*(2.d0*phi(lo(1)-1,lo(2),k,n) - phi(lo(1)-1,lo(2)+1,k,n))
                end if

                if (mask(hi(1)+1,lo(2)-1,k) .gt. 0) then
                   phi(hi(1)+1,lo(2)-1,k,n) = & ! Southeast
                        0.5d0*(2.d0*phi(hi(1),lo(2)-1,k,n) - phi(hi(1)-1,lo(2)-1,k,n)) + &
                        0.5d0*(2.d0*phi(hi(1)+1,lo(2),k,n) - phi(hi(1)+1,lo(2)+1,k,n))
                end if
             end do
          case (yhi_dir)
             do k = lo(3), hi(3)
                if (mask(lo(1)-1,hi(2)+1,k) .gt. 0) then
                   phi(lo(1)-1,hi(2)+1,k,n) = & ! Northwest
                        0.5d0*(2.d0*phi(lo(1),hi(2)+1,k,n) - phi(lo(1)+1,hi(2)+1,k,n)) + &
                        0.5d0*(2.d0*phi(lo(1)-1,hi(2),k,n) - phi(lo(1)-1,hi(2)-1,k,n))
                end if

                if (mask(hi(1)+1,hi(2)+1,k) .gt. 0) then
                   phi(hi(1)+1,hi(2)+1,k,n) = & ! Northeast
                        0.5d0*(2.d0*phi(hi(1),hi(2)+1,k,n) - phi(hi(1)-1,hi(2)+1,k,n)) + &
                        0.5d0*(2.d0*phi(hi(1)+1,hi(2),k,n) - phi(hi(1)+1,hi(2)-1,k,n))
                end if
             end do
          end select
       end if
#endif

    end do
  end subroutine amrex_mllinop_apply_bc


  subroutine amrex_mllinop_comp_interp_coef0 (lo, hi, &
       den, dlo, dhi, &
       mask, mlo, mhi, &
       cdir, bct, bcl, maxorder, dxinv, &
       nc) bind(c,name='amrex_mllinop_comp_interp_coef0')
    integer, dimension(3), intent(in) :: lo, hi, dlo, dhi, mlo, mhi
    integer, value, intent(in) :: cdir, bct, maxorder, nc
    real(amrex_real), value, intent(in) :: bcl
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) ::  den(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nc)
    integer         , intent(in   ) :: mask(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k,idim,lenx,m,n
    real(amrex_real) ::    x(-1:maxorder-2)
    real(amrex_real) :: coef(-1:maxorder-2)
    real(amrex_real), parameter :: xInt = -0.5D0
    real(amrex_real) :: c0
    
    do n = 1, nc
       if (bct == LO_NEUMANN) then

          select case (cdir)
          case (xlo_dir)
             do    k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   den(lo(1),j,k,n) = 1.d0
                end do
             end do
          case (xhi_dir)
             do    k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   den(hi(1),j,k,n) = 1.d0
                end do
             end do
#if (AMREX_SPACEDIM >= 2)
          case (ylo_dir)
             do    k = lo(3), hi(3)
                do i = lo(1), hi(1)
                   den(i,lo(2),k,n) = 1.d0
                end do
             end do
          case (yhi_dir)
             do    k = lo(3), hi(3)
                do i = lo(1), hi(1)
                   den(i,hi(2),k,n) = 1.d0
                end do
             end do
#if (AMREX_SPACEDIM == 3)
          case (zlo_dir)
             do    j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   den(i,j,lo(3),n) = 1.d0
                end do
             end do
          case (zhi_dir)
             do    j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   den(i,j,hi(3),n) = 1.d0
                end do
             end do
#endif
#endif
          end select

       else if (bct == LO_REFLECT_ODD .or. bct == LO_DIRICHLET) then

          if (bct == LO_REFLECT_ODD) then
             c0 = 1.d0
          else
             idim = mod(cdir,amrex_spacedim) + 1 ! cdir starts with 0; idim starts with 1
             lenx = MIN(hi(idim)-lo(idim), maxorder-2)

             x(-1) = -bcl*dxinv(idim)
             do m=0,maxorder-2
                x(m) = m + 0.5D0
             end do

             call polyInterpCoeff(xInt, x, lenx+2, coef)

             c0 = coef(0)
          end if

          select case (cdir)
          case (xlo_dir)
             do    k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   if (mask(lo(1)-1,j,k) .gt. 0) then
                      den(lo(1),j,k,n) = c0
                   else
                      den(lo(1),j,k,n) = 0.d0
                   end if
                end do
             end do
          case (xhi_dir)
             do    k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   if (mask(hi(1)+1,j,k) .gt. 0) then
                      den(hi(1),j,k,n) = c0
                   else
                      den(hi(1),j,k,n) = 0.d0
                   end if
                end do
             end do
#if (AMREX_SPACEDIM >= 2)
          case (ylo_dir)
             do    k = lo(3), hi(3)
                do i = lo(1), hi(1)
                   if (mask(i,lo(2)-1,k) .gt. 0) then
                      den(i,lo(2),k,n) = c0
                   else
                      den(i,lo(2),k,n) = 0.d0
                   end if
                end do
             end do
          case (yhi_dir)
             do    k = lo(3), hi(3)
                do i = lo(1), hi(1)
                   if (mask(i,hi(2)+1,k) .gt. 0) then
                      den(i,hi(2),k,n) = c0
                   else
                      den(i,hi(2),k,n) = 0.d0
                   end if
                end do
             end do
#if (AMREX_SPACEDIM == 3)
          case (zlo_dir)
             do    j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if (mask(i,j,lo(3)-1) .gt. 0) then
                      den(i,j,lo(3),n) = c0
                   else
                      den(i,j,lo(3),n) = 0.d0
                   end if
                end do
             end do
          case (zhi_dir)
             do    j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if (mask(i,j,hi(3)+1) .gt. 0) then
                      den(i,j,hi(3),n) = c0
                   else
                      den(i,j,hi(3),n) = 0.d0
                   end if
                end do
             end do
#endif
#endif
          end select

       end if
    end do
    
  end subroutine amrex_mllinop_comp_interp_coef0


  subroutine amrex_mllinop_apply_metric (lo, hi, &
       d, dlo, dhi, &
       r, rlo, rhi, &
       nc) &
       bind(c,name='amrex_mllinop_apply_metric')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), rlo, rhi
    integer, intent(in), value :: nc
    real(amrex_real), intent(inout) :: d(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nc)
    real(amrex_real), intent(in) :: r(rlo:rhi)

    integer :: i,j,k,n
    do n = 1, nc
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                d(i,j,k,n) = d(i,j,k,n) * r(i)
             end do
          end do
       end do
    end do
  end subroutine amrex_mllinop_apply_metric
  
end module amrex_mllinop_nd_module
