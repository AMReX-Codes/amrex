
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
  
end module warpx_module
