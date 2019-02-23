
module amrex_mllinop_nd_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real, amrex_spacedim
  use amrex_lo_util_module, only : polyInterpCoeff
  use amrex_lo_bctypes_module
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
  public :: amrex_mllinop_apply_bc

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
       if (bct == amrex_lo_neumann .or. bct == amrex_lo_reflect_odd) then

          if (bct == amrex_lo_neumann) then
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

       else if (bct == amrex_lo_dirichlet) then

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
#if (AMREX_SPACEDIM > 1)
       if (cross .eq. 0) then
          ! The iteration over faces is always in the order of xlo, ylo, zlo, xhi, yhi and zhi.

          ! Corners in XY plane
          if (cdir==xlo_dir .or. cdir==ylo_dir) then
             do k = lo(3), hi(3)
                phi(lo(1)-1,lo(2)-1,k,n) = & 
                     0.5d0*(2.d0*phi(lo(1),lo(2)-1,k,n) - phi(lo(1)+1,lo(2)-1,k,n)) + &
                     0.5d0*(2.d0*phi(lo(1)-1,lo(2),k,n) - phi(lo(1)-1,lo(2)+1,k,n))
             end do
          end if
          if (cdir==xhi_dir .or. cdir==ylo_dir) then
             do k = lo(3), hi(3)
                phi(hi(1)+1,lo(2)-1,k,n) = & 
                     0.5d0*(2.d0*phi(hi(1),lo(2)-1,k,n) - phi(hi(1)-1,lo(2)-1,k,n)) + &
                     0.5d0*(2.d0*phi(hi(1)+1,lo(2),k,n) - phi(hi(1)+1,lo(2)+1,k,n))
             end do
          end if
          if (cdir==xlo_dir .or. cdir==yhi_dir) then
             do k = lo(3), hi(3)
                phi(lo(1)-1,hi(2)+1,k,n) = & 
                     0.5d0*(2.d0*phi(lo(1),hi(2)+1,k,n) - phi(lo(1)+1,hi(2)+1,k,n)) + &
                     0.5d0*(2.d0*phi(lo(1)-1,hi(2),k,n) - phi(lo(1)-1,hi(2)-1,k,n))
             end do
          end if
          if (cdir==xhi_dir .or. cdir==yhi_dir) then
             do k = lo(3), hi(3)
                phi(hi(1)+1,hi(2)+1,k,n) = & 
                     0.5d0*(2.d0*phi(hi(1),hi(2)+1,k,n) - phi(hi(1)-1,hi(2)+1,k,n)) + &
                     0.5d0*(2.d0*phi(hi(1)+1,hi(2),k,n) - phi(hi(1)+1,hi(2)-1,k,n))
             end do
          end if
#if (AMREX_SPACEDIM > 2)
          ! Corners in YZ plane
          if (cdir==zlo_dir .or. cdir==ylo_dir) then
             do i = lo(1), hi(1)
                phi(i,lo(2)-1,lo(3)-1,n) = & 
                     0.5d0*(2.d0*phi(i,lo(2)-1,lo(3),n) - phi(i,lo(2)-1,lo(3)+1,n)) + &
                     0.5d0*(2.d0*phi(i,lo(2),lo(3)-1,n) - phi(i,lo(2)+1,lo(3)-1,n))
             end do
          end if
          if (cdir==zhi_dir .or. cdir==ylo_dir) then
             do i = lo(1), hi(1)
                phi(i,lo(2)-1,hi(3)+1,n) = & 
                     0.5d0*(2.d0*phi(i,lo(2)-1,hi(3),n) - phi(i,lo(2)-1,hi(3)-1,n)) + &
                     0.5d0*(2.d0*phi(i,lo(2),hi(3)+1,n) - phi(i,lo(2)+1,hi(3)+1,n))
             end do
          end if
          if (cdir==zlo_dir .or. cdir==yhi_dir) then
             do i = lo(1), hi(1)
                phi(i,hi(2)+1,lo(3)-1,n) = & 
                     0.5d0*(2.d0*phi(i,hi(2)+1,lo(3),n) - phi(i,hi(2)+1,lo(3)+1,n)) + &
                     0.5d0*(2.d0*phi(i,hi(2),lo(3)-1,n) - phi(i,hi(2)-1,lo(3)-1,n))
             end do
          end if
          if (cdir==zhi_dir .or. cdir==yhi_dir) then
             do i = lo(1), hi(1)
                phi(i,hi(2)+1,hi(3)+1,n) = & 
                     0.5d0*(2.d0*phi(i,hi(2)+1,hi(3),n) - phi(i,hi(2)+1,hi(3)-1,n)) + &
                     0.5d0*(2.d0*phi(i,hi(2),hi(3)+1,n) - phi(i,hi(2)-1,hi(3)+1,n))
             end do
          end if
          ! Corners in XZ plane
          if (cdir==xlo_dir .or. cdir==zlo_dir) then
             do j = lo(2), hi(2)
                phi(lo(1)-1,j,lo(3)-1,n) = & 
                     0.5d0*(2.d0*phi(lo(1),  j,lo(3)-1,n) - phi(lo(1)+1,j,lo(3)-1,n)) + &
                     0.5d0*(2.d0*phi(lo(1)-1,j,lo(3),n) - phi(lo(1)-1,j,lo(3)+1,n))
             end do
          end if
          if (cdir==xhi_dir .or. cdir==zlo_dir) then
             do j = lo(2), hi(2)
                phi(hi(1)+1,j,lo(3)-1,n) = & 
                     0.5d0*(2.d0*phi(hi(1),j,lo(3)-1,n) - phi(hi(1)-1,j,lo(3)-1,n)) + &
                     0.5d0*(2.d0*phi(hi(1)+1,j,lo(3),n) - phi(hi(1)+1,j,lo(3)+1,n))
             end do
          end if
          if (cdir==xlo_dir .or. cdir==zhi_dir) then
             do j = lo(2), hi(2)
                phi(lo(1)-1,j,hi(3)+1,n) = & 
                     0.5d0*(2.d0*phi(lo(1),j,hi(3)+1,n) - phi(lo(1)+1,j,hi(3)+1,n)) + &
                     0.5d0*(2.d0*phi(lo(1)-1,j,hi(3),n) - phi(lo(1)-1,j,hi(3)-1,n))
             end do
          end if
          if (cdir==xhi_dir .or. cdir==zhi_dir) then
             do j = lo(2), hi(2)
                phi(hi(1)+1,j,hi(3)+1,n) = & 
                     0.5d0*(2.d0*phi(hi(1),j,hi(3)+1,n) - phi(hi(1)-1,j,hi(3)+1,n)) + &
                     0.5d0*(2.d0*phi(hi(1)+1,j,hi(3),n) - phi(hi(1)+1,j,hi(3)-1,n))
             end do
          end if
#endif
       end if
#endif

    end do
  end subroutine amrex_mllinop_apply_bc
  
end module amrex_mllinop_nd_module
