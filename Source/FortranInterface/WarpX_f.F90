
module warpx_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real, amrex_spacedim

  implicit none

contains

  subroutine warpx_copy_attribs(np, xp, yp, zp, uxp, uyp, uzp, &
       xpold, ypold, zpold, uxpold, uypold, uzpold) &
       bind(c,name='warpx_copy_attribs')
    integer(c_long),  intent(in)    :: np
    real(amrex_real), intent(in)    :: xp(np)
    real(amrex_real), intent(in)    :: yp(np)
    real(amrex_real), intent(in)    :: zp(np)
    real(amrex_real), intent(in)    :: uxp(np)
    real(amrex_real), intent(in)    :: uyp(np)
    real(amrex_real), intent(in)    :: uzp(np)
    real(amrex_real), intent(inout) :: xpold(np)
    real(amrex_real), intent(inout) :: ypold(np)
    real(amrex_real), intent(inout) :: zpold(np)
    real(amrex_real), intent(inout) :: uxpold(np)
    real(amrex_real), intent(inout) :: uypold(np)
    real(amrex_real), intent(inout) :: uzpold(np)

    integer n

    do n = 1, np

       xpold(n) = xp(n)
       ypold(n) = yp(n)
       zpold(n) = zp(n)

       uxpold(n) = uxp(n)
       uypold(n) = uyp(n)
       uzpold(n) = uzp(n)

    end do

  end subroutine warpx_copy_attribs

  subroutine warpx_compute_E (lo, hi, &
       phi, phlo, phhi, &
       Ex,  Exlo, Exhi, &
       Ey,  Eylo, Eyhi, &
       Ez,  Ezlo, Ezhi, &
       dx) bind(c,name='warpx_compute_E')
    integer(c_int), intent(in) :: lo(3), hi(3), phlo(3), phhi(3), Exlo(3), Exhi(3),  &
         Eylo(3), Eyhi(3), Ezlo(3), Ezhi(3)
    real(amrex_real), intent(in)  :: dx(3)
    real(amrex_real), intent(in   ) :: phi(phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))
    real(amrex_real), intent(inout) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),Exlo(3):Exhi(3))
    real(amrex_real), intent(inout) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),Eylo(3):Eyhi(3))
    real(amrex_real), intent(inout) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),Ezlo(3):Ezhi(3))

    integer :: i, j, k
    real(amrex_real) :: dxinv(3)

    dxinv = 1.0 / dx

    do k    = lo(3), hi(3)
       do j = lo(2), hi(2)

          do i = lo(1), hi(1)-1
             Ex(i,j,k) = dxinv(1) * (phi(i,j,k) - phi(i+1,j,k))
          end do

          if (j < hi(2)) then
             do i = lo(1), hi(1)
                Ey(i,j,k) = dxinv(2) * (phi(i,j,k) - phi(i,j+1,k))
             end do
          end if

          if (k < hi(3)) then
             do i = lo(1), hi(1)
                Ez(i,j,k) = dxinv(3) * (phi(i,j,k) - phi(i,j,k+1))
             end do
          end if

       end do
    end do

  end subroutine warpx_compute_E

  subroutine warpx_compute_divb_3d (lo, hi, divB, dlo, dhi, &
       Bx, xlo, xhi, By, ylo, yhi, Bz, zlo, zhi, dx) &
       bind(c, name='warpx_compute_divb_3d')
    integer, intent(in) :: lo(3),hi(3),dlo(3),dhi(3),xlo(3),xhi(3),ylo(3),yhi(3),zlo(3),zhi(3)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(in   ) :: Bx  (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    real(amrex_real), intent(in   ) :: By  (ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3))
    real(amrex_real), intent(in   ) :: Bz  (zlo(1):zhi(1),zlo(2):zhi(2),zlo(3):zhi(3))
    real(amrex_real), intent(inout) :: divB(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))

    integer :: i,j,k
    real(amrex_real) :: dxinv(3)

    dxinv = 1.d0/dx

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             divB(i,j,k) = dxinv(1) * (Bx(i+1,j  ,k  ) - Bx(i,j,k)) &
                  +        dxinv(2) * (By(i  ,j+1,k  ) - By(i,j,k)) &
                  +        dxinv(3) * (Bz(i  ,j  ,k+1) - Bz(i,j,k))
          end do
       end do
    end do
  end subroutine warpx_compute_divb_3d


  subroutine warpx_compute_divb_2d (lo, hi, divB, dlo, dhi, &
       Bx, xlo, xhi, By, ylo, yhi, Bz, zlo, zhi, dx) &
       bind(c, name='warpx_compute_divb_2d')
    integer, intent(in) :: lo(2),hi(2),dlo(2),dhi(2),xlo(2),xhi(2),ylo(2),yhi(2),zlo(2),zhi(2)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(in   ) :: Bx  (xlo(1):xhi(1),xlo(2):xhi(2))
    real(amrex_real), intent(in   ) :: By  (ylo(1):yhi(1),ylo(2):yhi(2))
    real(amrex_real), intent(in   ) :: Bz  (zlo(1):zhi(1),zlo(2):zhi(2))
    real(amrex_real), intent(inout) :: divB(dlo(1):dhi(1),dlo(2):dhi(2))

    integer :: i,k
    real(amrex_real) :: dxinv(3)

    dxinv = 1.d0/dx

    do    k = lo(2), hi(2)
       do i = lo(1), hi(1)
          divB(i,k) = dxinv(1) * (Bx(i+1,k  ) - Bx(i,k)) &
               +      dxinv(3) * (Bz(i  ,k+1) - Bz(i,k))
       end do
    end do
  end subroutine warpx_compute_divb_2d


  subroutine warpx_compute_dive_3d (lo, hi, dive, dlo, dhi, &
       Ex, xlo, xhi, Ey, ylo, yhi, Ez, zlo, zhi, dx) &
       bind(c, name='warpx_compute_dive_3d')
    integer, intent(in) :: lo(3),hi(3),dlo(3),dhi(3),xlo(3),xhi(3),ylo(3),yhi(3),zlo(3),zhi(3)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(in   ) :: Ex  (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    real(amrex_real), intent(in   ) :: Ey  (ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3))
    real(amrex_real), intent(in   ) :: Ez  (zlo(1):zhi(1),zlo(2):zhi(2),zlo(3):zhi(3))
    real(amrex_real), intent(inout) :: dive(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))

    integer :: i,j,k
    real(amrex_real) :: dxinv(3)

    dxinv = 1.d0/dx

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dive(i,j,k) = dxinv(1) * (Ex(i,j,k) - Ex(i-1,j,k)) &
                  +        dxinv(2) * (Ey(i,j,k) - Ey(i,j-1,k)) &
                  +        dxinv(3) * (Ez(i,j,k) - Ez(i,j,k-1))
          end do
       end do
    end do
  end subroutine warpx_compute_dive_3d


  subroutine warpx_compute_dive_2d (lo, hi, dive, dlo, dhi, &
       Ex, xlo, xhi, Ey, ylo, yhi, Ez, zlo, zhi, dx) &
       bind(c, name='warpx_compute_dive_2d')
    integer, intent(in) :: lo(2),hi(2),dlo(2),dhi(2),xlo(2),xhi(2),ylo(2),yhi(2),zlo(2),zhi(2)
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(in   ) :: Ex  (xlo(1):xhi(1),xlo(2):xhi(2))
    real(amrex_real), intent(in   ) :: Ey  (ylo(1):yhi(1),ylo(2):yhi(2))
    real(amrex_real), intent(in   ) :: Ez  (zlo(1):zhi(1),zlo(2):zhi(2))
    real(amrex_real), intent(inout) :: dive(dlo(1):dhi(1),dlo(2):dhi(2))

    integer :: i,k
    real(amrex_real) :: dxinv(3)

    dxinv = 1.d0/dx

    do    k = lo(2), hi(2)
       do i = lo(1), hi(1)
          dive(i,k) = dxinv(1) * (Ex(i,k) - Ex(i-1,k)) &
               +      dxinv(3) * (Ez(i,k) - Ez(i,k-1))
       end do
    end do
  end subroutine warpx_compute_dive_2d


  subroutine warpx_compute_dive_rz (lo, hi, dive, dlo, dhi, &
       Ex, xlo, xhi, Ey, ylo, yhi, Ez, zlo, zhi, dx, rmin) &
       bind(c, name='warpx_compute_dive_rz')
    integer, intent(in) :: lo(2),hi(2),dlo(2),dhi(2),xlo(2),xhi(2),ylo(2),yhi(2),zlo(2),zhi(2)
    real(amrex_real), intent(in) :: dx(3), rmin
    real(amrex_real), intent(in   ) :: Ex  (xlo(1):xhi(1),xlo(2):xhi(2))
    real(amrex_real), intent(in   ) :: Ey  (ylo(1):yhi(1),ylo(2):yhi(2))
    real(amrex_real), intent(in   ) :: Ez  (zlo(1):zhi(1),zlo(2):zhi(2))
    real(amrex_real), intent(inout) :: dive(dlo(1):dhi(1),dlo(2):dhi(2))

    integer :: i,k
    real(amrex_real) :: dxinv(3)
    real(amrex_real) :: ru, rd

    dxinv = 1.d0/dx

    do    k = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (i == 0 .and. rmin == 0.) then
             ! the bulk equation diverges on axis
             ! (due to the 1/r terms). The following expressions regularize
             ! these divergences.
             dive(i,k) = 4.*dxinv(1) * Ex(i,k) &
                         + dxinv(3) * (Ez(i,k) - Ez(i,k-1))
          else
             ru = 1.d0 + 0.5d0/(rmin*dxinv(1) + i)
             rd = 1.d0 - 0.5d0/(rmin*dxinv(1) + i)
             dive(i,k) = dxinv(1) * (ru*Ex(i,k) - rd*Ex(i-1,k)) &
                       + dxinv(3) * (Ez(i,k) - Ez(i,k-1))
          end if
       end do
    end do
  end subroutine warpx_compute_dive_rz


  subroutine warpx_sync_current_2d (lo, hi, crse, clo, chi, fine, flo, fhi, dir) &
       bind(c, name='warpx_sync_current_2d')
    integer, intent(in) :: lo(2), hi(2), flo(2), fhi(2), clo(2), chi(2), dir
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1),flo(2):fhi(2))
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2))

    integer :: i,j,ii,jj

    if (dir == 0) then
       do j = lo(2), hi(2)
          jj = j*2
          do i = lo(1), hi(1)
             ii = i*2
             crse(i,j) = 0.25d0 * (fine(ii,jj) + fine(ii+1,jj) &
                  + 0.5d0*(fine(ii,jj-1) + fine(ii+1,jj-1) + fine(ii,jj+1) + fine(ii+1,jj+1)) )
          end do
       end do
    else if (dir == 2) then
       do j = lo(2), hi(2)
          jj = j*2
          do i = lo(1), hi(1)
             ii = i*2
             crse(i,j) = 0.25d0 * (fine(ii,jj) + fine(ii,jj+1) &
                  + 0.5d0*(fine(ii-1,jj) + fine(ii-1,jj+1) + fine(ii+1,jj) + fine(ii+1,jj+1)) )
          end do
       end do
    else
       do j = lo(2), hi(2)
          jj = j*2
          do i = lo(1), hi(1)
             ii = i*2
             crse(i,j) = 0.25d0 * &
                     ( fine(ii,jj) + 0.5d0 *(fine(ii-1,jj  )+fine(ii+1,jj  ) &
                     &                     + fine(ii  ,jj-1)+fine(ii  ,jj+1)) &
                     &             + 0.25d0*(fine(ii-1,jj-1)+fine(ii+1,jj-1) &
                     &                     + fine(ii-1,jj+1)+fine(ii+1,jj+1)) )
          end do
       end do
    end if
  end subroutine warpx_sync_current_2d

  subroutine warpx_sync_current_3d (lo, hi, crse, clo, chi, fine, flo, fhi, dir) &
       bind(c, name='warpx_sync_current_3d')
    integer, intent(in) :: lo(3), hi(3), flo(3), fhi(3), clo(3), chi(3), dir
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))

    integer :: i,j,k, ii,jj,kk

    if (dir == 0) then
       do k = lo(3), hi(3)
          kk = k*2
          do j = lo(2), hi(2)
             jj = j*2
             do i = lo(1), hi(1)
                ii = i*2
                crse(i,j,k) = 0.125d0* &
                     ( fine(ii  ,jj,kk) + 0.5d0 *(fine(ii  ,jj-1,kk  )+fine(ii  ,jj+1,kk  ) &
                     &                          + fine(ii  ,jj  ,kk-1)+fine(ii  ,jj  ,kk+1)) &
                     &                  + 0.25d0*(fine(ii  ,jj-1,kk-1)+fine(ii  ,jj+1,kk-1) &
                     &                          + fine(ii  ,jj-1,kk+1)+fine(ii  ,jj+1,kk+1)) &
                     + fine(ii+1,jj,kk) + 0.5d0 *(fine(ii+1,jj-1,kk  )+fine(ii+1,jj+1,kk  ) &
                     &                          + fine(ii+1,jj  ,kk-1)+fine(ii+1,jj  ,kk+1)) &
                     &                  + 0.25d0*(fine(ii+1,jj-1,kk-1)+fine(ii+1,jj+1,kk-1) &
                     &                          + fine(ii+1,jj-1,kk+1)+fine(ii+1,jj+1,kk+1)) )
             end do
          end do
       end do
    else if (dir == 1) then
       do k = lo(3), hi(3)
          kk = k*2
          do j = lo(2), hi(2)
             jj = j*2
             do i = lo(1), hi(1)
                ii = i*2
                crse(i,j,k) = 0.125d0* &
                     ( fine(ii,jj  ,kk) + 0.5d0 *(fine(ii-1,jj  ,kk  )+fine(ii+1,jj  ,kk  ) &
                     &                          + fine(ii  ,jj  ,kk-1)+fine(ii  ,jj  ,kk+1)) &
                     &                  + 0.25d0*(fine(ii-1,jj  ,kk-1)+fine(ii+1,jj  ,kk-1) &
                     &                          + fine(ii-1,jj  ,kk+1)+fine(ii+1,jj  ,kk+1)) &
                     + fine(ii,jj+1,kk) + 0.5d0 *(fine(ii-1,jj+1,kk  )+fine(ii+1,jj+1,kk  ) &
                     &                          + fine(ii  ,jj+1,kk-1)+fine(ii  ,jj+1,kk+1)) &
                     &                  + 0.25d0*(fine(ii-1,jj+1,kk-1)+fine(ii+1,jj+1,kk-1) &
                     &                          + fine(ii-1,jj+1,kk+1)+fine(ii+1,jj+1,kk+1)) )
             end do
          end do
       end do
    else
       do k = lo(3), hi(3)
          kk = k*2
          do j = lo(2), hi(2)
             jj = j*2
             do i = lo(1), hi(1)
                ii = i*2
                crse(i,j,k) = 0.125d0* &
                     ( fine(ii,jj,kk  ) + 0.5d0 *(fine(ii-1,jj  ,kk  )+fine(ii+1,jj  ,kk  ) &
                     &                          + fine(ii  ,jj-1,kk  )+fine(ii  ,jj+1,kk  )) &
                     &                  + 0.25d0*(fine(ii-1,jj-1,kk  )+fine(ii+1,jj-1,kk  ) &
                     &                          + fine(ii-1,jj+1,kk  )+fine(ii+1,jj+1,kk  )) &
                     + fine(ii,jj,kk+1) + 0.5d0 *(fine(ii-1,jj  ,kk+1)+fine(ii+1,jj  ,kk+1) &
                     &                          + fine(ii  ,jj-1,kk+1)+fine(ii  ,jj+1,kk+1)) &
                     &                  + 0.25d0*(fine(ii-1,jj-1,kk+1)+fine(ii+1,jj-1,kk+1) &
                     &                          + fine(ii-1,jj+1,kk+1)+fine(ii+1,jj+1,kk+1)) )
             end do
          end do
       end do
    end if
  end subroutine warpx_sync_current_3d


  subroutine warpx_sync_rho_2d (lo, hi, crse, clo, chi, fine, flo, fhi, nc) &
       bind(c, name='warpx_sync_rho_2d')
    integer, intent(in) :: lo(2), hi(2), flo(2), fhi(2), clo(2), chi(2), nc
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1),flo(2):fhi(2),nc)
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2),nc)

    integer :: i,j,ii,jj,m

    do m = 1, nc
       do j = lo(2), hi(2)
          jj = j*2
          do i = lo(1), hi(1)
             ii = i*2
             crse(i,j,m) = 0.25d0 * &
                  ( fine(ii,jj,m) + 0.5d0 *(fine(ii-1,jj  ,m)+fine(ii+1,jj  ,m) &
                  &                       + fine(ii  ,jj-1,m)+fine(ii  ,jj+1,m)) &
                  &               + 0.25d0*(fine(ii-1,jj-1,m)+fine(ii+1,jj-1,m) &
                  &                       + fine(ii-1,jj+1,m)+fine(ii+1,jj+1,m)) )
          end do
       end do
    end do
  end subroutine warpx_sync_rho_2d

  subroutine warpx_sync_rho_3d (lo, hi, crse, clo, chi, fine, flo, fhi, nc) &
       bind(c, name='warpx_sync_rho_3d')
    integer, intent(in) :: lo(3), hi(3), flo(3), fhi(3), clo(3), chi(3), nc
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3),nc)
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3),nc)

    integer :: i,j,k,ii,jj,kk,m

    do m = 1, nc
       do k = lo(3), hi(3)
          kk = k*2
          do j = lo(2), hi(2)
             jj = j*2
             do i = lo(1), hi(1)
                ii = i*2
                crse(i,j,k,m) = 0.125d0 * &
                     (fine(ii,jj,kk,m) + 0.5d0  *(fine(ii-1,jj  ,kk  ,m)+fine(ii+1,jj  ,kk  ,m) &
                     &                          + fine(ii  ,jj-1,kk  ,m)+fine(ii  ,jj+1,kk  ,m) &
                     &                          + fine(ii  ,jj  ,kk-1,m)+fine(ii  ,jj  ,kk+1,m)) &
                     &                 + 0.25d0 *(fine(ii-1,jj-1,kk  ,m)+fine(ii+1,jj-1,kk  ,m) &
                     &                          + fine(ii-1,jj+1,kk  ,m)+fine(ii+1,jj+1,kk  ,m) &
                     &                          + fine(ii-1,jj  ,kk-1,m)+fine(ii+1,jj  ,kk-1,m) &
                     &                          + fine(ii-1,jj  ,kk+1,m)+fine(ii+1,jj  ,kk+1,m) &
                     &                          + fine(ii  ,jj-1,kk-1,m)+fine(ii  ,jj+1,kk-1,m) &
                     &                          + fine(ii  ,jj-1,kk+1,m)+fine(ii  ,jj+1,kk+1,m)) &
                     &                 + 0.125d0*(fine(ii-1,jj-1,kk-1,m)+fine(ii-1,jj-1,kk+1,m) &
                     &                          + fine(ii-1,jj+1,kk-1,m)+fine(ii-1,jj+1,kk+1,m) &
                     &                          + fine(ii+1,jj-1,kk-1,m)+fine(ii+1,jj-1,kk+1,m) &
                     &                          + fine(ii+1,jj+1,kk-1,m)+fine(ii+1,jj+1,kk+1,m)))
             end do
          end do
       end do
    end do

  end subroutine warpx_sync_rho_3d

  subroutine warpx_build_buffer_masks (lo, hi, msk, mlo, mhi, gmsk, glo, ghi, ng) &
       bind(c, name='warpx_build_buffer_masks')
    integer, dimension(3), intent(in) :: lo, hi, mlo, mhi, glo, ghi
    integer, intent(in   ) :: gmsk(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3))
    integer, intent(inout) ::  msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    integer, intent(in) :: ng

    integer :: i,j,k

    if (amrex_spacedim .eq. 2) then

       k = lo(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (any(gmsk(i-ng:i+ng,j-ng:j+ng,k).eq.0)) then
                msk(i,j,k) = 0
             else
                msk(i,j,k) = 1
             end if
          end do
       end do

    else

       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (any(gmsk(i-ng:i+ng,j-ng:j+ng,k-ng:k+ng).eq.0)) then
                   msk(i,j,k) = 0
                else
                   msk(i,j,k) = 1
                end if
             end do
          end do
       end do

    end if

  end subroutine warpx_build_buffer_masks

end module warpx_module
