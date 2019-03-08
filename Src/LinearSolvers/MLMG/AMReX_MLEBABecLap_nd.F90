module amrex_mlebabeclap_nd_module
  
  use amrex_error_module
  use amrex_fort_module, only : amrex_real, amrex_spacedim
  use amrex_constants_module, only : zero, one, half
  use amrex_lo_util_module, only : polyInterpCoeff
  use amrex_lo_bctypes_module, only : amrex_lo_neumann, amrex_lo_dirichlet, amrex_lo_reflect_odd
  use amrex_ebcellflag_module, only: is_regular_cell, is_covered_cell
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
  public :: amrex_mlebabeclap_apply_bc, amrex_eb_copy_dirichlet

contains

  subroutine amrex_mlebabeclap_apply_bc (lo, hi, phi, hlo, hhi, flag, flo, fhi, &
       apx, xlo, xhi, &
#if (AMREX_SPACEDIM >= 2)
       apy, ylo, yhi, &
#endif
#if (AMREX_SPACEDIM == 3)
       apz, zlo, zhi, &
#endif
       mask, mlo, mhi, cdir, bct, bcl, bcval, blo, bhi, maxorder_in, dxinv, inhomog, nc) &
       bind(c,name='amrex_mlebabeclap_apply_bc')
    integer, dimension(3), intent(in) :: lo, hi, hlo, hhi, flo, fhi, mlo, mhi, blo, bhi, xlo, xhi
#if (AMREX_SPACEDIM >= 2)
    integer, dimension(3), intent(in) :: ylo, yhi
#endif
#if (AMREX_SPACEDIM == 3)
    integer, dimension(3), intent(in) :: zlo, zhi
#endif    
    integer, value, intent(in) :: cdir, bct, maxorder_in, inhomog, nc
    real(amrex_real), value, intent(in) :: bcl
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) ::  phi (hlo(1):hhi(1),hlo(2):hhi(2),hlo(3):hhi(3),nc)
    integer         , intent(in   ) :: flag (flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    integer         , intent(in   ) :: mask (mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    real(amrex_real), intent(in   ) :: bcval(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),nc)
    real(amrex_real), intent(in   ) :: apx  (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
#if (AMREX_SPACEDIM >= 2)
    real(amrex_real), intent(in   ) :: apy  (ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3))
#endif
#if (AMREX_SPACEDIM == 3)
    real(amrex_real), intent(in   ) :: apz  (zlo(1):zhi(1),zlo(2):zhi(2),zlo(3):zhi(3))
#endif
    
    integer :: i, j, k, idim, m, n, iorder, maxorder, glo(3), ghi(3)
    logical :: inhomogeneous
    real(amrex_real) ::    x(-1:2)
    real(amrex_real) :: coef(-1:2,2:4), coef2(-2:1,2:4) ! maxorder can be 2, 3, or 4
    real(amrex_real), parameter :: xInt = -half
    real(amrex_real) :: fac

    glo = lo
    ghi = hi
    do i = 1, amrex_spacedim
       glo(i) = glo(i)-1
       ghi(i) = ghi(i)+1
    end do
    
    ! note that the mask in this subroutine is different from masks in bndry registers
    ! 1 means valid data, 0 means invalid data
    
    inhomogeneous = (inhomog .ne. 0)

    do n = 1, nc
       if (bct == amrex_lo_neumann .or. bct == amrex_lo_reflect_odd) then

          if (bct == amrex_lo_neumann) then
             fac = one
          else
             fac = -one
          end if

          select case (cdir)
          case (xlo_dir)
             do    k = glo(3), ghi(3)
                do j = glo(2), ghi(2)
                   if (mask(lo(1)-1,j,k) .eq. 0 .and. mask(lo(1),j,k) .eq. 1) then
                      phi(lo(1)-1,j,k,n) = fac*phi(lo(1),j,k,n)
                   end if
                end do
             end do
          case (xhi_dir)
             do    k = glo(3), ghi(3)
                do j = glo(2), ghi(2)
                   if (mask(hi(1)+1,j,k) .eq. 0 .and. mask(hi(1),j,k) .eq. 1) then
                      phi(hi(1)+1,j,k,n) = fac*phi(hi(1),j,k,n)
                   end if
                end do
             end do
#if (AMREX_SPACEDIM >= 2)
          case (ylo_dir)  
             do    k = glo(3), ghi(3)
                do i = glo(1), ghi(1)
                   if (mask(i,lo(2)-1,k) .eq. 0 .and. mask(i,lo(2),k) .eq. 1) then
                      phi(i,lo(2)-1,k,n) = fac*phi(i,lo(2),k,n)
                   end if
                end do
             end do
          case (yhi_dir)
             do    k = glo(3), ghi(3)
                do i = glo(1), ghi(1)
                   if (mask(i,hi(2)+1,k) .eq. 0 .and. mask(i,hi(2),k) .eq. 1) then
                      phi(i,hi(2)+1,k,n) = fac*phi(i,hi(2),k,n)
                   end if
                end do
             end do
#if (AMREX_SPACEDIM == 3)
          case (zlo_dir)
             do    j = glo(2), ghi(2)
                do i = glo(1), ghi(1)
                   if (mask(i,j,lo(3)-1) .eq. 0 .and. mask(i,j,lo(3)) .eq. 1) then
                      phi(i,j,lo(3)-1,n) = fac*phi(i,j,lo(3),n)
                   end if
                end do
             end do
          case (zhi_dir)
             do    j = glo(2), ghi(2)
                do i = glo(1), ghi(1)
                   if (mask(i,j,hi(3)+1) .eq. 0 .and. mask(i,j,hi(3)) .eq. 1) then
                      phi(i,j,hi(3)+1,n) = fac*phi(i,j,hi(3),n)
                   end if
                end do
             end do
#endif
#endif
          end select

       else if (bct == amrex_lo_dirichlet) then

          idim = mod(cdir,amrex_spacedim) + 1 ! cdir starts with 0; idim starts with 1
          maxorder = min(maxorder_in, hi(idim)-lo(idim)+2)

          x(-1) = -bcl*dxinv(idim)
          do m=0,2
             x(m) = m + 0.5D0
          end do

          do iorder = 2, maxorder
             call polyInterpCoeff(xInt, x, iorder, coef(:,iorder))
             do m = -(iorder-2), 1
                coef2(m,iorder) = coef(-m,iorder)
             end do
          end do

          select case (cdir)
          case (xlo_dir)
             do    k = glo(3), ghi(3)
                do j = glo(2), ghi(2)
                   if (mask(lo(1)-1,j,k) .eq. 0 .and. mask(lo(1),j,k) .eq. 1) then
                      iorder = 1
                      do m = 0, maxorder-2
                         if (apx(lo(1)+m,j,k) .gt. zero) then
                            iorder = iorder + 1
                         else
                            exit
                         end if
                      end do
                      if (iorder .eq. 1) then
                         if (inhomogeneous) then
                            phi(lo(1)-1,j,k,n) = bcval(lo(1)-1,j,k,n)
                         else
                            phi(lo(1)-1,j,k,n) = zero
                         end if
                      else
                         phi(lo(1)-1,j,k,n) = sum(phi(lo(1):lo(1)+(iorder-2),j,k,n) &
                              * coef(0:(iorder-2),iorder))
                         if (inhomogeneous) then
                            phi(lo(1)-1,j,k,n) = phi(lo(1)-1,j,k,n) &
                                 + bcval(lo(1)-1,j,k,n)*coef(-1,iorder)
                         end if
                      end if
                   end if
                end do
             end do
          case (xhi_dir)
             do    k = glo(3), ghi(3)
                do j = glo(2), ghi(2)
                   if (mask(hi(1)+1,j,k) .eq. 0 .and. mask(hi(1),j,k) .eq. 1) then
                      iorder = 1
                      do m = 0, maxorder-2
                         if (apx(hi(1)+1-m,j,k) .gt. zero) then
                            iorder = iorder + 1
                         else
                            exit
                         end if
                      end do
                      if (iorder .eq. 1) then
                         if (inhomogeneous) then
                            phi(hi(1)+1,j,k,n) = bcval(hi(1)+1,j,k,n)
                         else
                            phi(hi(1)+1,j,k,n) = zero
                         end if
                      else
                         phi(hi(1)+1,j,k,n) = sum(phi(hi(1)-(iorder-2):hi(1),j,k,n) &
                              * coef2(-(iorder-2):0,iorder))
                         if (inhomogeneous) then
                            phi(hi(1)+1,j,k,n) = phi(hi(1)+1,j,k,n) &
                                 + bcval(hi(1)+1,j,k,n)*coef2(1,iorder)
                         end if
                      end if
                   end if
                end do
             end do
#if (AMREX_SPACEDIM >= 2)
          case (ylo_dir)
             do    k = glo(3), ghi(3)
                do i = glo(1), ghi(1)
                   if (mask(i,lo(2)-1,k) .eq. 0 .and. mask(i,lo(2),k) .eq. 1) then
                      iorder = 1
                      do m = 0, maxorder-2
                         if (apy(i,lo(2)+m,k) .gt. zero) then
                            iorder = iorder + 1
                         else
                            exit
                         end if
                      end do
                      if (iorder .eq. 1) then
                         if (inhomogeneous) then
                            phi(i,lo(2)-1,k,n) = bcval(i,lo(2)-1,k,n)
                         else
                            phi(i,lo(2)-1,k,n) = zero
                         end if
                      else
                         phi(i,lo(2)-1,k,n) = sum(phi(i,lo(2):lo(2)+(iorder-2),k,n) &
                              * coef(0:(iorder-2),iorder))
                         if (inhomogeneous) then
                            phi(i,lo(2)-1,k,n) = phi(i,lo(2)-1,k,n) &
                                 + bcval(i,lo(2)-1,k,n)*coef(-1,iorder)
                         end if
                      end if
                   end if
                end do
             end do
          case (yhi_dir)
             do    k = glo(3), ghi(3)
                do i = glo(1), ghi(1)
                   if (mask(i,hi(2)+1,k) .eq. 0 .and. mask(i,hi(2),k) .eq. 1) then
                      iorder = 1
                      do m = 0, maxorder-2
                         if (apy(i,hi(2)+1-m,k) .gt. zero) then
                            iorder = iorder + 1
                         else
                            exit
                         end if
                      end do
                      if (iorder .eq. 1) then
                         if (inhomogeneous) then
                            phi(i,hi(2)+1,k,n) = bcval(i,hi(2)+1,k,n)
                         else
                            phi(i,hi(2)+1,k,n) = zero
                         end if
                      else
                         phi(i,hi(2)+1,k,n) = sum(phi(i,hi(2)-(iorder-2):hi(2),k,n) &
                              * coef2(-(iorder-2):0,iorder))
                         if (inhomogeneous) then
                            phi(i,hi(2)+1,k,n) = phi(i,hi(2)+1,k,n) &
                                 + bcval(i,hi(2)+1,k,n)*coef2(1,iorder)
                         end if
                      end if
                   end if
                end do
             end do
#if (AMREX_SPACEDIM == 3)
          case (zlo_dir)
             do    j = glo(2), ghi(2)
                do i = glo(1), ghi(1)
                   if (mask(i,j,lo(3)-1) .eq. 0 .and. mask(i,j,lo(3)) .eq. 1) then
                      iorder = 1
                      do m = 0, maxorder-2
                         if (apz(i,j,lo(3)+m) .gt. zero) then
                            iorder = iorder +1
                         else
                            exit
                         end if
                      end do
                      if (iorder .eq. 1) then
                         if (inhomogeneous) then
                            phi(i,j,lo(3)-1,n) = bcval(i,j,lo(3)-1,n)
                         else
                            phi(i,j,lo(3)-1,n) = zero
                         end if
                      else
                         phi(i,j,lo(3)-1,n) = sum(phi(i,j,lo(3):lo(3)+(iorder-2),n)&
                              * coef(0:(iorder-2),iorder))
                         if (inhomogeneous) then
                            phi(i,j,lo(3)-1,n) = phi(i,j,lo(3)-1,n) &
                                 + bcval(i,j,lo(3)-1,n)*coef(-1,iorder)
                         end if
                      end if
                   end if
                end do
             end do
          case (zhi_dir)
             do    j = glo(2), ghi(2)
                do i = glo(1), ghi(1)
                   if (mask(i,j,hi(3)+1) .eq. 0 .and. mask(i,j,hi(3)) .eq. 1) then
                      iorder = 1
                      do m = 0, maxorder-2
                         if (apz(i,j,hi(3)+1-m) .gt. zero) then
                            iorder = iorder + 1
                         else
                            exit
                         end if
                      end do
                      if (iorder .eq. 1) then
                         if (inhomogeneous) then
                            phi(i,j,hi(3)+1,n) = bcval(i,j,hi(3)+1,n)
                         else
                            phi(i,j,hi(3)+1,n) = zero
                         end if
                      else
                         phi(i,j,hi(3)+1,n) = sum(phi(i,j,hi(3)-(iorder-2):hi(3),n) &
                              * coef2(-(iorder-2):0,iorder))
                         if (inhomogeneous) then
                            phi(i,j,hi(3)+1,n) = phi(i,j,hi(3)+1,n) &
                                 + bcval(i,j,hi(3)+1,n)*coef2(1,iorder)
                         end if
                      end if
                   end if
                end do
             end do
#endif
#endif
          end select

       else
          call amrex_error("amrex_mlebabeclap_nd_module: unknown bc");
       end if
       
    end do

  end subroutine amrex_mlebabeclap_apply_bc

  subroutine amrex_eb_copy_dirichlet (lo, hi, phi, phlo, phhi, phiin, pilo, pihi, beta, blo, bhi, &
       betain, bilo, bihi, flag, flo, fhi) bind(c,name='amrex_eb_copy_dirichlet')

    use amrex_ebcellflag_module, only : is_single_valued_cell
    integer, dimension(3), intent(in) :: lo, hi, phlo, phhi, pilo, pihi, blo, bhi, bilo, bihi, flo, fhi
    real(amrex_real), intent(inout) :: phi   (phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))
    real(amrex_real), intent(in   ) :: phiin (pilo(1):pihi(1),pilo(2):pihi(2),pilo(3):pihi(3))
    real(amrex_real), intent(inout) :: beta  ( blo(1): bhi(1), blo(2): bhi(2), blo(3): bhi(3))
    real(amrex_real), intent(in   ) :: betain(bilo(1):bihi(1),bilo(2):bihi(2),bilo(3):bihi(3))
    integer         , intent(in   ) :: flag  ( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3))

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (is_single_valued_cell(flag(i,j,k))) then
                phi(i,j,k) = phiin(i,j,k)
                beta(i,j,k) = betain(i,j,k)
             else
                phi(i,j,k) = zero
                beta(i,j,k) = zero
             end if
          end do
       end do
    end do
  end subroutine amrex_eb_copy_dirichlet

  subroutine amrex_eb_homog_dirichlet (lo, hi, phi, phlo, phhi,  beta, blo, bhi, &
       betain, bilo, bihi, flag, flo, fhi) bind(c,name='amrex_eb_homog_dirichlet')

    use amrex_ebcellflag_module, only : is_single_valued_cell
    integer, dimension(3), intent(in) :: lo, hi, phlo, phhi, blo, bhi, bilo, bihi, flo, fhi
    real(amrex_real), intent(inout) :: phi   (phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))
    real(amrex_real), intent(inout) :: beta  ( blo(1): bhi(1), blo(2): bhi(2), blo(3): bhi(3))
    real(amrex_real), intent(in   ) :: betain(bilo(1):bihi(1),bilo(2):bihi(2),bilo(3):bihi(3))
    integer         , intent(in   ) :: flag  ( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3))

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (is_single_valued_cell(flag(i,j,k))) then
                beta(i,j,k) = betain(i,j,k)
             else
                beta(i,j,k) = zero
             end if

             phi(i,j,k)  = zero

          end do
       end do
    end do
  end subroutine amrex_eb_homog_dirichlet
  
end module amrex_mlebabeclap_nd_module
