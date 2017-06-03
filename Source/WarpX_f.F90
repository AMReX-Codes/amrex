
module warpx_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real

  implicit none

contains

  subroutine warpx_compute_E (lo, hi, &
       phi, phlo, phhi, &
       Ex,  Exlo, Exhi, &
       Ey,  Eylo, Eyhi, &
       Ez,  Ezlo, Ezhi, &
       dx) bind(c,name='warpx_compute_E')
    integer(c_int), intent(in) :: lo(3), hi(3), phlo(3), phhi(3), Exlo(3), Exhi(3),  &
         Eylo(3), Eyhi(3), Ezlo(3), Ezhi(3)
    real(amrex_real), intent(in) :: dx(3)
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
       Bx, xlo, xhi, Bz, zlo, zhi, dx) &
       bind(c, name='warpx_compute_divb_2d')
    integer, intent(in) :: lo(2),hi(2),dlo(2),dhi(2),xlo(2),xhi(2),zlo(2),zhi(2)
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(in   ) :: Bx  (xlo(1):xhi(1),xlo(2):xhi(2))
    real(amrex_real), intent(in   ) :: Bz  (zlo(1):zhi(1),zlo(2):zhi(2))
    real(amrex_real), intent(inout) :: divB(dlo(1):dhi(1),dlo(2):dhi(2))

    integer :: i,k
    real(amrex_real) :: dxinv(2)

    dxinv = 1.d0/dx

    do    k = lo(2), hi(2)
       do i = lo(1), hi(1)
          divB(i,k) = dxinv(1) * (Bx(i+1,k  ) - Bx(i,k)) &
               +      dxinv(2) * (Bz(i  ,k+1) - Bz(i,k))
       end do
    end do
  end subroutine warpx_compute_divb_2d


  subroutine warpx_push_pml_bvec_3d (xlo, xhi, ylo, yhi, zlo, zhi, &
       &                             Ex, Exlo, Exhi, &
       &                             Ey, Eylo, Eyhi, &
       &                             Ez, Ezlo, Ezhi, &
       &                             Bx, Bxlo, Bxhi, &
       &                             By, Bylo, Byhi, &
       &                             Bz, Bzlo, Bzhi, &
       &                             sigx1, sigx1_lo, sigx1_hi, &
       &                             sigx2, sigx2_lo, sigx2_hi, &
       &                             sigy1, sigy1_lo, sigy1_hi, &
       &                             sigy2, sigy2_lo, sigy2_hi, &
       &                             sigz1, sigz1_lo, sigz1_hi, &
       &                             sigz2, sigz2_lo, sigz2_hi) &
       bind(c,name='warpx_push_pml_bvec_3d')
    integer, intent(in) :: xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
         Exlo(3), Exhi(3), Eylo(3), Eyhi(3), Ezlo(3), Ezhi(3), &
         Bxlo(3), Bxhi(3), Bylo(3), Byhi(3), Bzlo(3), Bzhi(3)
    integer, intent(in), value :: sigx1_lo, sigx1_hi, sigx2_lo, sigx2_hi
    integer, intent(in), value :: sigy1_lo, sigy1_hi, sigy2_lo, sigy2_hi
    integer, intent(in), value :: sigz1_lo, sigz1_hi, sigz2_lo, sigz2_hi
    real(amrex_real), intent(in   ) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),Exlo(3):Exhi(3),2)
    real(amrex_real), intent(in   ) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),Eylo(3):Eyhi(3),2)
    real(amrex_real), intent(in   ) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),Ezlo(3):Ezhi(3),2)
    real(amrex_real), intent(inout) :: Bx (Bxlo(1):Bxhi(1),Bxlo(2):Bxhi(2),Bxlo(3):Bxhi(3),2)
    real(amrex_real), intent(inout) :: By (Bylo(1):Byhi(1),Bylo(2):Byhi(2),Bylo(3):Byhi(3),2)
    real(amrex_real), intent(inout) :: Bz (Bzlo(1):Bzhi(1),Bzlo(2):Bzhi(2),Bzlo(3):Bzhi(3),2)
    real(amrex_real), intent(in) :: sigx1(sigx1_lo:sigx1_hi)
    real(amrex_real), intent(in) :: sigx2(sigx2_lo:sigx2_hi)
    real(amrex_real), intent(in) :: sigy1(sigy1_lo:sigy1_hi)
    real(amrex_real), intent(in) :: sigy2(sigy2_lo:sigy2_hi)
    real(amrex_real), intent(in) :: sigz1(sigz1_lo:sigz1_hi)
    real(amrex_real), intent(in) :: sigz2(sigz2_lo:sigz2_hi)

    integer :: i, j, k

    do       k = xlo(3), xhi(3)
       do    j = xlo(2), xhi(2)
          do i = xlo(1), xhi(1)
             Bx(i,j,k,1) = sigy1(j)*Bx(i,j,k,1) - sigy2(j)*(Ez(i,j+1,k  ,1)+Ez(i,j+1,k  ,2) &
                  &                                        -Ez(i,j  ,k  ,1)-Ez(i,j  ,k  ,2))
             Bx(i,j,k,2) = sigz1(k)*Bx(i,j,k,2) + sigz2(k)*(Ey(i,j  ,k+1,1)+Ey(i,j  ,k+1,2) &
                  &                                        -Ey(i,j  ,k  ,1)-Ey(i,j  ,k  ,2))
          end do
       end do
    end do

    do       k = ylo(3), yhi(3)
       do    j = ylo(2), yhi(2)
          do i = ylo(1), yhi(1)
             By(i,j,k,1) = sigz1(k)*By(i,j,k,1) - sigz2(k)*(Ex(i  ,j,k+1,1)+Ex(i  ,j,k+1,2) &
                  &                                        -Ex(i  ,j,k  ,1)-Ex(i  ,j,k  ,2))
             By(i,j,k,2) = sigx1(i)*By(i,j,k,2) + sigx2(i)*(Ez(i+1,j,k  ,1)+Ez(i+1,j,k  ,2) &
                  &                                        -Ez(i  ,j,k  ,1)-Ez(i  ,j,k  ,2))
          end do
       end do
    end do

    do       k = zlo(3), zhi(3)
       do    j = zlo(2), zhi(2)
          do i = zlo(1), zhi(1)
             Bz(i,j,k,1) = sigx1(i)*Bz(i,j,k,1) - sigx2(i)*(Ey(i+1,j  ,k,1)+Ey(i+1,j  ,k,2) &
                  &                                        -Ey(i  ,j  ,k,1)-Ey(i  ,j  ,k,2))
             Bz(i,j,k,2) = sigy1(j)*Bz(i,j,k,2) + sigy2(j)*(Ex(i  ,j+1,k,1)+Ex(i  ,j+1,k,2) &
                  &                                        -Ex(i  ,j  ,k,1)-Ex(i  ,j  ,k,2))
          end do
       end do
    end do

  end subroutine warpx_push_pml_bvec_3d

  subroutine warpx_push_pml_evec_3d (xlo, xhi, ylo, yhi, zlo, zhi, &
       &                             Ex, Exlo, Exhi, &
       &                             Ey, Eylo, Eyhi, &
       &                             Ez, Ezlo, Ezhi, &
       &                             Bx, Bxlo, Bxhi, &
       &                             By, Bylo, Byhi, &
       &                             Bz, Bzlo, Bzhi, &
       &                             sigx1, sigx1_lo, sigx1_hi, &
       &                             sigx2, sigx2_lo, sigx2_hi, &
       &                             sigy1, sigy1_lo, sigy1_hi, &
       &                             sigy2, sigy2_lo, sigy2_hi, &
       &                             sigz1, sigz1_lo, sigz1_hi, &
       &                             sigz2, sigz2_lo, sigz2_hi) &
       bind(c,name='warpx_push_pml_evec_3d')
    integer, intent(in) :: xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
         Exlo(3), Exhi(3), Eylo(3), Eyhi(3), Ezlo(3), Ezhi(3), &
         Bxlo(3), Bxhi(3), Bylo(3), Byhi(3), Bzlo(3), Bzhi(3)
    integer, intent(in), value :: sigx1_lo, sigx1_hi, sigx2_lo, sigx2_hi
    integer, intent(in), value :: sigy1_lo, sigy1_hi, sigy2_lo, sigy2_hi
    integer, intent(in), value :: sigz1_lo, sigz1_hi, sigz2_lo, sigz2_hi
    real(amrex_real), intent(inout) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),Exlo(3):Exhi(3),2)
    real(amrex_real), intent(inout) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),Eylo(3):Eyhi(3),2)
    real(amrex_real), intent(inout) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),Ezlo(3):Ezhi(3),2)
    real(amrex_real), intent(in   ) :: Bx (Bxlo(1):Bxhi(1),Bxlo(2):Bxhi(2),Bxlo(3):Bxhi(3),2)
    real(amrex_real), intent(in   ) :: By (Bylo(1):Byhi(1),Bylo(2):Byhi(2),Bylo(3):Byhi(3),2)
    real(amrex_real), intent(in   ) :: Bz (Bzlo(1):Bzhi(1),Bzlo(2):Bzhi(2),Bzlo(3):Bzhi(3),2)
    real(amrex_real), intent(in) :: sigx1(sigx1_lo:sigx1_hi)
    real(amrex_real), intent(in) :: sigx2(sigx2_lo:sigx2_hi)
    real(amrex_real), intent(in) :: sigy1(sigy1_lo:sigy1_hi)
    real(amrex_real), intent(in) :: sigy2(sigy2_lo:sigy2_hi)
    real(amrex_real), intent(in) :: sigz1(sigz1_lo:sigz1_hi)
    real(amrex_real), intent(in) :: sigz2(sigz2_lo:sigz2_hi)

    integer :: i, j, k

    do       k = xlo(3), xhi(3)
       do    j = xlo(2), xhi(2)
          do i = xlo(1), xhi(1)
             Ex(i,j,k,1) = sigy1(j)*Ex(i,j,k,1) + sigy2(j)*(Bz(i,j  ,k  ,1)+Bz(i,j  ,k  ,2) &
                  &                                        -Bz(i,j-1,k  ,1)-Bz(i,j-1,k  ,2))
             Ex(i,j,k,2) = sigz1(k)*Ex(i,j,k,2) - sigz2(k)*(By(i,j  ,k  ,1)+By(i,j  ,k  ,2) &
                  &                                        -By(i,j  ,k-1,1)-By(i,j  ,k-1,2))
          end do
       end do
    end do

    do       k = ylo(3), yhi(3)
       do    j = ylo(2), yhi(2)
          do i = ylo(1), yhi(1)
             Ey(i,j,k,1) = sigz1(k)*Ey(i,j,k,1) + sigz2(k)*(Bx(i  ,j,k  ,1)+Bx(i  ,j,k  ,2) &
                  &                                        -Bx(i  ,j,k-1,1)-Bx(i  ,j,k-1,2))
             Ey(i,j,k,2) = sigx1(i)*Ey(i,j,k,2) - sigx2(i)*(Bz(i  ,j,k  ,1)+Bz(i  ,j,k  ,2) &
                  &                                        -Bz(i-1,j,k  ,1)-Bz(i-1,j,k  ,2))
          end do
       end do
    end do

    do       k = zlo(3), zhi(3)
       do    j = zlo(2), zhi(2)
          do i = zlo(1), zhi(1)
             Ez(i,j,k,1) = sigx1(i)*Ez(i,j,k,1) + sigx2(i)*(By(i  ,j  ,k,1)+By(i  ,j  ,k,2) &
                  &                                        -By(i-1,j  ,k,1)-By(i-1,j  ,k,2))
             Ez(i,j,k,2) = sigy1(j)*Ez(i,j,k,2) - sigy2(j)*(Bx(i  ,j  ,k,1)+Bx(i  ,j  ,k,2) &
                  &                                        -Bx(i  ,j-1,k,1)-Bx(i  ,j-1,k,2))
          end do
       end do
    end do

  end subroutine warpx_push_pml_evec_3d

  subroutine warpx_push_pml_bvec_2d (xlo, xhi, ylo, yhi, zlo, zhi, &
       &                             Ex, Exlo, Exhi, &
       &                             Ey, Eylo, Eyhi, &
       &                             Ez, Ezlo, Ezhi, &
       &                             Bx, Bxlo, Bxhi, &
       &                             By, Bylo, Byhi, &
       &                             Bz, Bzlo, Bzhi, &
       &                             sigx1, sigx1_lo, sigx1_hi, &
       &                             sigx2, sigx2_lo, sigx2_hi, &
       &                             sigz1, sigz1_lo, sigz1_hi, &
       &                             sigz2, sigz2_lo, sigz2_hi) &
       bind(c,name='warpx_push_pml_bvec_2d')
    integer, intent(in) :: xlo(2), xhi(2), ylo(2), yhi(2), zlo(2), zhi(2), &
         Exlo(2), Exhi(2), Eylo(2), Eyhi(2), Ezlo(2), Ezhi(2), &
         Bxlo(2), Bxhi(2), Bylo(2), Byhi(2), Bzlo(2), Bzhi(2)
    integer, intent(in), value :: sigx1_lo, sigx1_hi, sigx2_lo, sigx2_hi
    integer, intent(in), value :: sigz1_lo, sigz1_hi, sigz2_lo, sigz2_hi
    real(amrex_real), intent(in   ) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),2)
    real(amrex_real), intent(in   ) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),2)
    real(amrex_real), intent(in   ) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),2)
    real(amrex_real), intent(inout) :: Bx (Bxlo(1):Bxhi(1),Bxlo(2):Bxhi(2),2)
    real(amrex_real), intent(inout) :: By (Bylo(1):Byhi(1),Bylo(2):Byhi(2),2)
    real(amrex_real), intent(inout) :: Bz (Bzlo(1):Bzhi(1),Bzlo(2):Bzhi(2),2)
    real(amrex_real), intent(in) :: sigx1(sigx1_lo:sigx1_hi)
    real(amrex_real), intent(in) :: sigx2(sigx2_lo:sigx2_hi)
    real(amrex_real), intent(in) :: sigz1(sigz1_lo:sigz1_hi)
    real(amrex_real), intent(in) :: sigz2(sigz2_lo:sigz2_hi)

    integer :: i, k

    do    k = xlo(2), xhi(2)
       do i = xlo(1), xhi(1)
          Bx(i,k,2) = sigz1(k)*Bx(i,k,2) + sigz2(k)*(Ey(i,k+1,1)+Ey(i,k+1,2) &
               &                                    -Ey(i,k  ,1)-Ey(i,k  ,2))
       end do
    end do

    do    k = ylo(2), yhi(2)
       do i = ylo(1), yhi(1)
          By(i,k,1) = sigz1(k)*By(i,k,1) - sigz2(k)*(Ex(i  ,k+1,1)+Ex(i  ,k+1,2) &
               &                                    -Ex(i  ,k  ,1)-Ex(i  ,k  ,2))
          By(i,k,2) = sigx1(i)*By(i,k,2) + sigx2(i)*(Ez(i+1,k  ,1)+Ez(i+1,k  ,2) &
               &                                    -Ez(i  ,k  ,1)-Ez(i  ,k  ,2))
       end do
    end do

    do    k = zlo(2), zhi(2)
       do i = zlo(1), zhi(1)
          Bz(i,k,1) = sigx1(i)*Bz(i,k,1) - sigx2(i)*(Ey(i+1,k,1)+Ey(i+1,k,2) &
               &                                    -Ey(i  ,k,1)-Ey(i  ,k,2))
       end do
    end do

  end subroutine warpx_push_pml_bvec_2d

  subroutine warpx_push_pml_evec_2d (xlo, xhi, ylo, yhi, zlo, zhi, &
       &                             Ex, Exlo, Exhi, &
       &                             Ey, Eylo, Eyhi, &
       &                             Ez, Ezlo, Ezhi, &
       &                             Bx, Bxlo, Bxhi, &
       &                             By, Bylo, Byhi, &
       &                             Bz, Bzlo, Bzhi, &
       &                             sigx1, sigx1_lo, sigx1_hi, &
       &                             sigx2, sigx2_lo, sigx2_hi, &
       &                             sigz1, sigz1_lo, sigz1_hi, &
       &                             sigz2, sigz2_lo, sigz2_hi) &
       bind(c,name='warpx_push_pml_evec_2d')
    integer, intent(in) :: xlo(2), xhi(2), ylo(2), yhi(2), zlo(2), zhi(2), &
         Exlo(2), Exhi(2), Eylo(2), Eyhi(2), Ezlo(2), Ezhi(2), &
         Bxlo(2), Bxhi(2), Bylo(2), Byhi(2), Bzlo(2), Bzhi(2)
    integer, intent(in), value :: sigx1_lo, sigx1_hi, sigx2_lo, sigx2_hi
    integer, intent(in), value :: sigz1_lo, sigz1_hi, sigz2_lo, sigz2_hi
    real(amrex_real), intent(inout) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),2)
    real(amrex_real), intent(inout) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),2)
    real(amrex_real), intent(inout) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),2)
    real(amrex_real), intent(in   ) :: Bx (Bxlo(1):Bxhi(1),Bxlo(2):Bxhi(2),2)
    real(amrex_real), intent(in   ) :: By (Bylo(1):Byhi(1),Bylo(2):Byhi(2),2)
    real(amrex_real), intent(in   ) :: Bz (Bzlo(1):Bzhi(1),Bzlo(2):Bzhi(2),2)
    real(amrex_real), intent(in) :: sigx1(sigx1_lo:sigx1_hi)
    real(amrex_real), intent(in) :: sigx2(sigx2_lo:sigx2_hi)
    real(amrex_real), intent(in) :: sigz1(sigz1_lo:sigz1_hi)
    real(amrex_real), intent(in) :: sigz2(sigz2_lo:sigz2_hi)

    integer :: i, k

    do    k = xlo(2), xhi(2)
       do i = xlo(1), xhi(1)
          Ex(i,k,2) = sigz1(k)*Ex(i,k,2) - sigz2(k)*(By(i,k  ,1)+By(i,k  ,2) &
               &                                    -By(i,k-1,1)-By(i,k-1,2))
       end do
    end do

    do    k = ylo(2), yhi(2)
       do i = ylo(1), yhi(1)
          Ey(i,k,1) = sigz1(k)*Ey(i,k,1) + sigz2(k)*(Bx(i  ,k  ,1)+Bx(i  ,k  ,2) &
               &                                    -Bx(i  ,k-1,1)-Bx(i  ,k-1,2))
          Ey(i,k,2) = sigx1(i)*Ey(i,k,2) - sigx2(i)*(Bz(i  ,k  ,1)+Bz(i  ,k  ,2) &
               &                                    -Bz(i-1,k  ,1)-Bz(i-1,k  ,2))
       end do
    end do

    do    k = zlo(2), zhi(2)
       do i = zlo(1), zhi(1)
          Ez(i,k,1) = sigx1(i)*Ez(i,k,1) + sigx2(i)*(By(i  ,k,1)+By(i  ,k,2) &
               &                                    -By(i-1,k,1)-By(i-1,k,2))
       end do
    end do

  end subroutine warpx_push_pml_evec_2d


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

end module warpx_module
