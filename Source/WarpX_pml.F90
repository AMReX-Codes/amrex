module warpx_pml_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real

  implicit none

contains

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
       &                             sigz2, sigz2_lo, sigz2_hi, &
       &                             dtsdx, dtsdy, dtsdz, solver_type) &
       bind(c,name='warpx_push_pml_bvec_3d')
    integer, intent(in) :: xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
         Exlo(3), Exhi(3), Eylo(3), Eyhi(3), Ezlo(3), Ezhi(3), &
         Bxlo(3), Bxhi(3), Bylo(3), Byhi(3), Bzlo(3), Bzhi(3), solver_type
    integer, intent(in), value :: sigx1_lo, sigx1_hi, sigx2_lo, sigx2_hi
    integer, intent(in), value :: sigy1_lo, sigy1_hi, sigy2_lo, sigy2_hi
    integer, intent(in), value :: sigz1_lo, sigz1_hi, sigz2_lo, sigz2_hi
    real(amrex_real), intent(in   ) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),Exlo(3):Exhi(3),3)
    real(amrex_real), intent(in   ) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),Eylo(3):Eyhi(3),3)
    real(amrex_real), intent(in   ) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),Ezlo(3):Ezhi(3),3)
    real(amrex_real), intent(inout) :: Bx (Bxlo(1):Bxhi(1),Bxlo(2):Bxhi(2),Bxlo(3):Bxhi(3),2)
    real(amrex_real), intent(inout) :: By (Bylo(1):Byhi(1),Bylo(2):Byhi(2),Bylo(3):Byhi(3),2)
    real(amrex_real), intent(inout) :: Bz (Bzlo(1):Bzhi(1),Bzlo(2):Bzhi(2),Bzlo(3):Bzhi(3),2)
    real(amrex_real), intent(in) :: sigx1(sigx1_lo:sigx1_hi)
    real(amrex_real), intent(in) :: sigx2(sigx2_lo:sigx2_hi)
    real(amrex_real), intent(in) :: sigy1(sigy1_lo:sigy1_hi)
    real(amrex_real), intent(in) :: sigy2(sigy2_lo:sigy2_hi)
    real(amrex_real), intent(in) :: sigz1(sigz1_lo:sigz1_hi)
    real(amrex_real), intent(in) :: sigz2(sigz2_lo:sigz2_hi)
    real(amrex_real), intent(in) :: dtsdx, dtsdy, dtsdz

    integer :: i, j, k
    
    real(amrex_real) :: delta, rx, ry, rz, betaxz, betaxy, betayx, betayz, betazx, betazy
    real(amrex_real) :: beta, alphax, alphay, alphaz, gammax, gammay, gammaz
  
    ! solver_type: 0=Yee; 1=CKC

    if (solver_type==0) then

        ! Yee push

        do       k = xlo(3), xhi(3)
           do    j = xlo(2), xhi(2)
              do i = xlo(1), xhi(1)
                 Bx(i,j,k,1) = sigy1(j)*Bx(i,j,k,1) - sigy2(j)*(Ez(i,j+1,k  ,1)+Ez(i,j+1,k  ,2)+Ez(i,j+1,k  ,3) &
                      &                                        -Ez(i,j  ,k  ,1)-Ez(i,j  ,k  ,2)-Ez(i,j  ,k  ,3))
                 Bx(i,j,k,2) = sigz1(k)*Bx(i,j,k,2) + sigz2(k)*(Ey(i,j  ,k+1,1)+Ey(i,j  ,k+1,2)+Ey(i,j  ,k+1,3) &
                      &                                        -Ey(i,j  ,k  ,1)-Ey(i,j  ,k  ,2)-Ey(i,j  ,k  ,3))
              end do
           end do
        end do

        do       k = ylo(3), yhi(3)
           do    j = ylo(2), yhi(2)
              do i = ylo(1), yhi(1)
                 By(i,j,k,1) = sigz1(k)*By(i,j,k,1) - sigz2(k)*(Ex(i  ,j,k+1,1)+Ex(i  ,j,k+1,2)+Ex(i  ,j,k+1,3) &
                      &                                        -Ex(i  ,j,k  ,1)-Ex(i  ,j,k  ,2)-Ex(i  ,j,k  ,3))
                 By(i,j,k,2) = sigx1(i)*By(i,j,k,2) + sigx2(i)*(Ez(i+1,j,k  ,1)+Ez(i+1,j,k  ,2)+Ez(i+1,j,k  ,3) &
                      &                                        -Ez(i  ,j,k  ,1)-Ez(i  ,j,k  ,2)-Ez(i  ,j,k  ,3))
              end do
           end do
        end do

        do       k = zlo(3), zhi(3)
           do    j = zlo(2), zhi(2)
              do i = zlo(1), zhi(1)
                 Bz(i,j,k,1) = sigx1(i)*Bz(i,j,k,1) - sigx2(i)*(Ey(i+1,j  ,k,1)+Ey(i+1,j  ,k,2)+Ey(i+1,j  ,k,3) &
                      &                                        -Ey(i  ,j  ,k,1)-Ey(i  ,j  ,k,2)-Ey(i  ,j  ,k,3))
                 Bz(i,j,k,2) = sigy1(j)*Bz(i,j,k,2) + sigy2(j)*(Ex(i  ,j+1,k,1)+Ex(i  ,j+1,k,2)+Ex(i  ,j+1,k,3) &
                      &                                        -Ex(i  ,j  ,k,1)-Ex(i  ,j  ,k,2)-Ex(i  ,j  ,k,3))
              end do
           end do
        end do
        
    else

        ! CKC push  

        ! computes coefficients according to Cowan - PRST-AB 16, 041303 (2013)
        delta = max(dtsdx,dtsdy,dtsdz)
        rx = (dtsdx/delta)**2
        ry = (dtsdy/delta)**2
        rz = (dtsdz/delta)**2
        beta = 0.125*(1.-rx*ry*rz/(ry*rz+rz*rx+rx*ry))
        betaxy = ry*beta
        betaxz = rz*beta
        betayx = rx*beta
        betayz = rz*beta
        betazx = rx*beta
        betazy = ry*beta
        gammax = ry*rz*(1./16.-0.125*ry*rz/(ry*rz+rz*rx+rx*ry))
        gammay = rx*rz*(1./16.-0.125*rx*rz/(ry*rz+rz*rx+rx*ry))
        gammaz = rx*ry*(1./16.-0.125*rx*ry/(ry*rz+rz*rx+rx*ry))
        alphax = 1. - 2.*betaxy - 2.* betaxz - 4.*gammax
        alphay = 1. - 2.*betayx - 2.* betayz - 4.*gammay
        alphaz = 1. - 2.*betazx - 2.* betazy - 4.*gammaz
    
        do       k = xlo(3), xhi(3)
           do    j = xlo(2), xhi(2)
              do i = xlo(1), xhi(1)
                 Bx(i,j,k,1) = sigy1(j)*Bx(i,j,k,1) - sigy2(j)*(alphay*(Ez(i  ,j+1,k  ,1)+Ez(i  ,j+1,k  ,2)+Ez(i  ,j+1,k  ,3)  &
                                                                       -Ez(i  ,j  ,k  ,1)-Ez(i  ,j  ,k  ,2)-Ez(i  ,j  ,k  ,3)) &
                                                               +betayx*(Ez(i+1,j+1,k  ,1)+Ez(i+1,j+1,k  ,2)+Ez(i+1,j+1,k  ,3)  &
                                                                       -Ez(i+1,j  ,k  ,1)-Ez(i+1,j  ,k  ,2)-Ez(i+1,j  ,k  ,3)  &
                                                                       +Ez(i-1,j+1,k  ,1)+Ez(i-1,j+1,k  ,2)+Ez(i-1,j+1,k  ,3)  &
                                                                       -Ez(i-1,j  ,k  ,1)-Ez(i-1,j  ,k  ,2)-Ez(i-1,j  ,k  ,3)) &
                                                               +betayz*(Ez(i  ,j+1,k+1,1)+Ez(i  ,j+1,k+1,2)+Ez(i  ,j+1,k+1,3)  &
                                                                       -Ez(i  ,j  ,k+1,1)-Ez(i  ,j  ,k+1,2)-Ez(i  ,j  ,k+1,3)  &
                                                                       +Ez(i  ,j+1,k-1,1)+Ez(i  ,j+1,k-1,2)+Ez(i  ,j+1,k-1,3)  &
                                                                       -Ez(i  ,j  ,k-1,1)-Ez(i  ,j  ,k-1,2)-Ez(i  ,j  ,k-1,3)) &
                                                               +gammay*(Ez(i+1,j+1,k+1,1)+Ez(i+1,j+1,k+1,2)+Ez(i+1,j+1,k+1,3)  &
                                                                       -Ez(i+1,j  ,k+1,1)-Ez(i+1,j  ,k+1,2)-Ez(i+1,j  ,k+1,3)  &
                                                                       +Ez(i-1,j+1,k+1,1)+Ez(i-1,j+1,k+1,2)+Ez(i-1,j+1,k+1,3)  &
                                                                       -Ez(i-1,j  ,k+1,1)-Ez(i-1,j  ,k+1,2)-Ez(i-1,j  ,k+1,3)  &
                                                                       +Ez(i+1,j+1,k-1,1)+Ez(i+1,j+1,k-1,2)+Ez(i+1,j+1,k-1,3)  &
                                                                       -Ez(i+1,j  ,k-1,1)-Ez(i+1,j  ,k-1,2)-Ez(i+1,j  ,k-1,3)  &
                                                                       +Ez(i-1,j+1,k-1,1)+Ez(i-1,j+1,k-1,2)+Ez(i-1,j+1,k-1,3)  &
                                                                       -Ez(i-1,j  ,k-1,1)-Ez(i-1,j  ,k-1,2)-Ez(i-1,j  ,k-1,3)))


                 Bx(i,j,k,2) = sigz1(k)*Bx(i,j,k,2) + sigz2(k)*(alphaz*(Ey(i  ,j  ,k+1,1)+Ey(i  ,j  ,k+1,2)+Ey(i  ,j  ,k+1,3)  &
                                                                       -Ey(i  ,j  ,k  ,1)-Ey(i  ,j  ,k  ,2)-Ey(i  ,j  ,k  ,3)) &
                                                               +betazx*(Ey(i+1,j  ,k+1,1)+Ey(i+1,j  ,k+1,2)+Ey(i+1,j  ,k+1,3)  &
                                                                       -Ey(i+1,j  ,k  ,1)-Ey(i+1,j  ,k  ,2)-Ey(i+1,j  ,k  ,3)  &
                                                                       +Ey(i-1,j  ,k+1,1)+Ey(i-1,j  ,k+1,2)+Ey(i-1,j  ,k+1,3)  &
                                                                       -Ey(i-1,j  ,k  ,1)-Ey(i-1,j  ,k  ,2)-Ey(i-1,j  ,k  ,3)) &
                                                               +betazy*(Ey(i  ,j+1,k+1,1)+Ey(i  ,j+1,k+1,2)+Ey(i  ,j+1,k+1,3)  &
                                                                       -Ey(i  ,j+1,k  ,1)-Ey(i  ,j+1,k  ,2)-Ey(i  ,j+1,k  ,3)  &
                                                                       +Ey(i  ,j-1,k+1,1)+Ey(i  ,j-1,k+1,2)+Ey(i  ,j-1,k+1,3)  &
                                                                       -Ey(i  ,j-1,k  ,1)-Ey(i  ,j-1,k  ,2)-Ey(i  ,j-1,k  ,3)) &
                                                               +gammaz*(Ey(i+1,j+1,k+1,1)+Ey(i+1,j+1,k+1,2)+Ey(i+1,j+1,k+1,3)  &
                                                                       -Ey(i+1,j+1,k  ,1)-Ey(i+1,j+1,k  ,2)-Ey(i+1,j+1,k  ,3)  &
                                                                       +Ey(i-1,j+1,k+1,1)+Ey(i-1,j+1,k+1,2)+Ey(i-1,j+1,k+1,3)  &
                                                                       -Ey(i-1,j+1,k  ,1)-Ey(i-1,j+1,k  ,2)-Ey(i-1,j+1,k  ,3)  &
                                                                       +Ey(i+1,j-1,k+1,1)+Ey(i+1,j-1,k+1,2)+Ey(i+1,j-1,k+1,3)  &
                                                                       -Ey(i+1,j-1,k  ,1)-Ey(i+1,j-1,k  ,2)-Ey(i+1,j-1,k  ,3)  &
                                                                       +Ey(i-1,j-1,k+1,1)+Ey(i-1,j-1,k+1,2)+Ey(i-1,j-1,k+1,3)  &
                                                                       -Ey(i-1,j-1,k  ,1)-Ey(i-1,j-1,k  ,2)-Ey(i-1,j-1,k  ,3)))
              end do
           end do
        end do

        do       k = ylo(3), yhi(3)
           do    j = ylo(2), yhi(2)
              do i = ylo(1), yhi(1)
                 By(i,j,k,1) = sigz1(k)*By(i,j,k,1) - sigz2(k)*(alphaz*(Ex(i  ,j  ,k+1,1)+Ex(i  ,j  ,k+1,2)+Ex(i  ,j  ,k+1,3)  &
                                                                       -Ex(i  ,j  ,k  ,1)-Ex(i  ,j  ,k  ,2)-Ex(i  ,j  ,k  ,3)) &
                                                               +betazx*(Ex(i+1,j  ,k+1,1)+Ex(i+1,j  ,k+1,2)+Ex(i+1,j  ,k+1,3)  &
                                                                       -Ex(i+1,j  ,k  ,1)-Ex(i+1,j  ,k  ,2)-Ex(i+1,j  ,k  ,3)  &
                                                                       +Ex(i-1,j  ,k+1,1)+Ex(i-1,j  ,k+1,2)+Ex(i-1,j  ,k+1,3)  &
                                                                       -Ex(i-1,j  ,k  ,1)-Ex(i-1,j  ,k  ,2)-Ex(i-1,j  ,k  ,3)) &
                                                               +betaxy*(Ex(i  ,j+1,k+1,1)+Ex(i  ,j+1,k+1,2)+Ex(i  ,j+1,k+1,3)  &
                                                                       -Ex(i  ,j+1,k  ,1)-Ex(i  ,j+1,k  ,2)-Ex(i  ,j+1,k  ,3)  &
                                                                       +Ex(i  ,j-1,k+1,1)+Ex(i  ,j-1,k+1,2)+Ex(i  ,j-1,k+1,3)  &
                                                                       -Ex(i  ,j-1,k  ,1)-Ex(i  ,j-1,k  ,2)-Ex(i  ,j-1,k  ,3)) &
                                                               +gammaz*(Ex(i+1,j+1,k+1,1)+Ex(i+1,j+1,k+1,2)+Ex(i+1,j+1,k+1,3)  &
                                                                       -Ex(i+1,j+1,k  ,1)-Ex(i+1,j+1,k  ,2)-Ex(i+1,j+1,k  ,3)  &
                                                                       +Ex(i-1,j+1,k+1,1)+Ex(i-1,j+1,k+1,2)+Ex(i-1,j+1,k+1,3)  &
                                                                       -Ex(i-1,j+1,k  ,1)-Ex(i-1,j+1,k  ,2)-Ex(i-1,j+1,k  ,3)  &
                                                                       +Ex(i+1,j-1,k+1,1)+Ex(i+1,j-1,k+1,2)+Ex(i+1,j-1,k+1,3)  &
                                                                       -Ex(i+1,j-1,k  ,1)-Ex(i+1,j-1,k  ,2)-Ex(i+1,j-1,k  ,3)  &
                                                                       +Ex(i-1,j-1,k+1,1)+Ex(i-1,j-1,k+1,2)+Ex(i-1,j-1,k+1,3)  &
                                                                       -Ex(i-1,j-1,k  ,1)-Ex(i-1,j-1,k  ,2)-Ex(i-1,j-1,k  ,3)))


                 By(i,j,k,2) = sigx1(i)*By(i,j,k,2) + sigx2(i)*(alphax*(Ez(i+1,j  ,k  ,1)+Ez(i+1,j  ,k  ,2)+Ez(i+1,j  ,k  ,3)  &
                                                                       -Ez(i  ,j  ,k  ,1)-Ez(i  ,j  ,k  ,2)-Ez(i  ,j  ,k  ,3)) &
                                                               +betaxy*(Ez(i+1,j+1,k  ,1)+Ez(i+1,j+1,k  ,2)+Ez(i+1,j+1,k  ,3)  &
                                                                       -Ez(i  ,j+1,k  ,1)-Ez(i  ,j+1,k  ,2)-Ez(i  ,j+1,k  ,3)  &
                                                                       +Ez(i+1,j-1,k  ,1)+Ez(i+1,j-1,k  ,2)+Ez(i+1,j-1,k  ,3)  &
                                                                       -Ez(i  ,j-1,k  ,1)-Ez(i  ,j-1,k  ,2)-Ez(i  ,j-1,k  ,3)) &
                                                               +betaxz*(Ez(i+1,j  ,k+1,1)+Ez(i+1,j  ,k+1,2)+Ez(i+1,j  ,k+1,3)  &
                                                                       -Ez(i  ,j  ,k+1,1)-Ez(i  ,j  ,k+1,2)-Ez(i  ,j  ,k+1,3)  &
                                                                       +Ez(i+1,j  ,k-1,1)+Ez(i+1,j  ,k-1,2)+Ez(i+1,j  ,k-1,3)  &
                                                                       -Ez(i  ,j  ,k-1,1)-Ez(i  ,j  ,k-1,2)-Ez(i  ,j  ,k-1,3)) &
                                                               +gammax*(Ez(i+1,j+1,k+1,1)+Ez(i+1,j+1,k+1,2)+Ez(i+1,j+1,k+1,3)  &
                                                                       -Ez(i  ,j+1,k+1,1)-Ez(i  ,j+1,k+1,2)-Ez(i  ,j+1,k+1,3)  &
                                                                       +Ez(i+1,j-1,k+1,1)+Ez(i+1,j-1,k+1,2)+Ez(i+1,j-1,k+1,3)  &
                                                                       -Ez(i  ,j-1,k+1,1)-Ez(i  ,j-1,k+1,2)-Ez(i  ,j-1,k+1,3)  &
                                                                       +Ez(i+1,j+1,k-1,1)+Ez(i+1,j+1,k-1,2)+Ez(i+1,j+1,k-1,3)  &
                                                                       -Ez(i  ,j+1,k-1,1)-Ez(i  ,j+1,k-1,2)-Ez(i  ,j+1,k-1,3)  &
                                                                       +Ez(i+1,j-1,k-1,1)+Ez(i+1,j-1,k-1,2)+Ez(i+1,j-1,k-1,3)  &
                                                                       -Ez(i  ,j-1,k-1,1)-Ez(i  ,j-1,k-1,2)-Ez(i  ,j-1,k-1,3)))
              end do
           end do
        end do

        do       k = zlo(3), zhi(3)
           do    j = zlo(2), zhi(2)
              do i = zlo(1), zhi(1)
                 Bz(i,j,k,1) = sigx1(i)*Bz(i,j,k,1) - sigx2(i)*(alphax*(Ey(i+1,j  ,k  ,1)+Ey(i+1,j  ,k  ,2)+Ey(i+1,j  ,k  ,3)  &
                                                                       -Ey(i  ,j  ,k  ,1)-Ey(i  ,j  ,k  ,2)-Ey(i  ,j  ,k  ,3)) &
                                                               +betaxy*(Ey(i+1,j+1,k  ,1)+Ey(i+1,j+1,k  ,2)+Ey(i+1,j+1,k  ,3)  &
                                                                       -Ey(i  ,j+1,k  ,1)-Ey(i  ,j+1,k  ,2)-Ey(i  ,j+1,k  ,3)  &
                                                                       +Ey(i+1,j-1,k  ,1)+Ey(i+1,j-1,k  ,2)+Ey(i+1,j-1,k  ,3)  &
                                                                       -Ey(i  ,j-1,k  ,1)-Ey(i  ,j-1,k  ,2)-Ey(i  ,j-1,k  ,3)) &
                                                               +betaxz*(Ey(i+1,j  ,k+1,1)+Ey(i+1,j  ,k+1,2)+Ey(i+1,j  ,k+1,3)  &
                                                                       -Ey(i  ,j  ,k+1,1)-Ey(i  ,j  ,k+1,2)-Ey(i  ,j  ,k+1,3)  &
                                                                       +Ey(i+1,j  ,k-1,1)+Ey(i+1,j  ,k-1,2)+Ey(i+1,j  ,k-1,3)  &
                                                                       -Ey(i  ,j  ,k-1,1)-Ey(i  ,j  ,k-1,2)-Ey(i  ,j  ,k-1,3)) &
                                                               +gammax*(Ey(i+1,j+1,k+1,1)+Ey(i+1,j+1,k+1,2)+Ey(i+1,j+1,k+1,3)  &
                                                                       -Ey(i  ,j+1,k+1,1)-Ey(i  ,j+1,k+1,2)-Ey(i  ,j+1,k+1,3)  &
                                                                       +Ey(i+1,j-1,k+1,1)+Ey(i+1,j-1,k+1,2)+Ey(i+1,j-1,k+1,3)  &
                                                                       -Ey(i  ,j-1,k+1,1)-Ey(i  ,j-1,k+1,2)-Ey(i  ,j-1,k+1,3)  &
                                                                       +Ey(i+1,j+1,k-1,1)+Ey(i+1,j+1,k-1,2)+Ey(i+1,j+1,k-1,3)  &
                                                                       -Ey(i  ,j+1,k-1,1)-Ey(i  ,j+1,k-1,2)-Ey(i  ,j+1,k-1,3)  &
                                                                       +Ey(i+1,j-1,k-1,1)+Ey(i+1,j-1,k-1,2)+Ey(i+1,j-1,k-1,3)  &
                                                                       -Ey(i  ,j-1,k-1,1)-Ey(i  ,j-1,k-1,2)-Ey(i  ,j-1,k-1,3)))

                 Bz(i,j,k,2) = sigy1(j)*Bz(i,j,k,2) + sigy2(j)*(alphay*(Ex(i  ,j+1,k  ,1)+Ex(i  ,j+1,k  ,2)+Ex(i  ,j+1,k  ,3)  &
                                                                       -Ex(i  ,j  ,k  ,1)-Ex(i  ,j  ,k  ,2)-Ex(i  ,j  ,k  ,3)) &
                                                               +betayx*(Ex(i+1,j+1,k  ,1)+Ex(i+1,j+1,k  ,2)+Ex(i+1,j+1,k  ,3)  &
                                                                       -Ex(i+1,j  ,k  ,1)-Ex(i+1,j  ,k  ,2)-Ex(i+1,j  ,k  ,3)  &
                                                                       +Ex(i-1,j+1,k  ,1)+Ex(i-1,j+1,k  ,2)+Ex(i-1,j+1,k  ,3)  &
                                                                       -Ex(i-1,j  ,k  ,1)-Ex(i-1,j  ,k  ,2)-Ex(i-1,j  ,k  ,3)) &
                                                               +betayz*(Ex(i  ,j+1,k+1,1)+Ex(i  ,j+1,k+1,2)+Ex(i  ,j+1,k+1,3)  &
                                                                       -Ex(i  ,j  ,k+1,1)-Ex(i  ,j  ,k+1,2)-Ex(i  ,j  ,k+1,3)  &
                                                                       +Ex(i  ,j+1,k-1,1)+Ex(i  ,j+1,k-1,2)+Ex(i  ,j+1,k-1,3)  &
                                                                       -Ex(i  ,j  ,k-1,1)-Ex(i  ,j  ,k-1,2)-Ex(i  ,j  ,k-1,3)) &
                                                               +gammay*(Ex(i+1,j+1,k+1,1)+Ex(i+1,j+1,k+1,2)+Ex(i+1,j+1,k+1,3)  &
                                                                       -Ex(i+1,j  ,k+1,1)-Ex(i+1,j  ,k+1,2)-Ex(i+1,j  ,k+1,3)  &
                                                                       +Ex(i-1,j+1,k+1,1)+Ex(i-1,j+1,k+1,2)+Ex(i-1,j+1,k+1,3)  &
                                                                       -Ex(i-1,j  ,k+1,1)-Ex(i-1,j  ,k+1,2)-Ex(i-1,j  ,k+1,3)  &
                                                                       +Ex(i+1,j+1,k-1,1)+Ex(i+1,j+1,k-1,2)+Ex(i+1,j+1,k-1,3)  &
                                                                       -Ex(i+1,j  ,k-1,1)-Ex(i+1,j  ,k-1,2)-Ex(i+1,j  ,k-1,3)  &
                                                                       +Ex(i-1,j+1,k-1,1)+Ex(i-1,j+1,k-1,2)+Ex(i-1,j+1,k-1,3)  &
                                                                       -Ex(i-1,j  ,k-1,1)-Ex(i-1,j  ,k-1,2)-Ex(i-1,j  ,k-1,3)))
              end do
           end do
        end do
        
    end if

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
       &                             sigz2, sigz2_lo, sigz2_hi, &
       &                             dtsdx, dtsdy, dtsdz, solver_type) &
       bind(c,name='warpx_push_pml_bvec_2d')
    integer, intent(in) :: xlo(2), xhi(2), ylo(2), yhi(2), zlo(2), zhi(2), &
         Exlo(2), Exhi(2), Eylo(2), Eyhi(2), Ezlo(2), Ezhi(2), &
         Bxlo(2), Bxhi(2), Bylo(2), Byhi(2), Bzlo(2), Bzhi(2), solver_type
    integer, intent(in), value :: sigx1_lo, sigx1_hi, sigx2_lo, sigx2_hi
    integer, intent(in), value :: sigz1_lo, sigz1_hi, sigz2_lo, sigz2_hi
    real(amrex_real), intent(in   ) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),3)
    real(amrex_real), intent(in   ) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),3)
    real(amrex_real), intent(in   ) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),3)
    real(amrex_real), intent(inout) :: Bx (Bxlo(1):Bxhi(1),Bxlo(2):Bxhi(2),2)
    real(amrex_real), intent(inout) :: By (Bylo(1):Byhi(1),Bylo(2):Byhi(2),2)
    real(amrex_real), intent(inout) :: Bz (Bzlo(1):Bzhi(1),Bzlo(2):Bzhi(2),2)
    real(amrex_real), intent(in) :: sigx1(sigx1_lo:sigx1_hi)
    real(amrex_real), intent(in) :: sigx2(sigx2_lo:sigx2_hi)
    real(amrex_real), intent(in) :: sigz1(sigz1_lo:sigz1_hi)
    real(amrex_real), intent(in) :: sigz2(sigz2_lo:sigz2_hi)
    real(amrex_real), intent(in) :: dtsdx, dtsdy, dtsdz

    integer :: i, k
    
    real(amrex_real) :: delta, rx, rz, betaxz, betazx, alphax, alphaz
  
    
    ! solver_type: 0=Yee; 1=CKC

    if (solver_type==0) then
  
    ! Yee push

        do    k = xlo(2), xhi(2)
           do i = xlo(1), xhi(1)
              Bx(i,k,2) = sigz1(k)*Bx(i,k,2) + sigz2(k)*(Ey(i,k+1,1)+Ey(i,k+1,2)+Ey(i,k+1,3) &
                   &                                    -Ey(i,k  ,1)-Ey(i,k  ,2)-Ey(i,k  ,3))
           end do
        end do

        do    k = ylo(2), yhi(2)
           do i = ylo(1), yhi(1)
              By(i,k,1) = sigz1(k)*By(i,k,1) - sigz2(k)*(Ex(i  ,k+1,1)+Ex(i  ,k+1,2)+Ex(i  ,k+1,3) &
                   &                                    -Ex(i  ,k  ,1)-Ex(i  ,k  ,2)-Ex(i  ,k  ,3))
              By(i,k,2) = sigx1(i)*By(i,k,2) + sigx2(i)*(Ez(i+1,k  ,1)+Ez(i+1,k  ,2)+Ez(i+1,k  ,3) &
                   &                                    -Ez(i  ,k  ,1)-Ez(i  ,k  ,2)-Ez(i  ,k  ,3))
           end do
        end do

        do    k = zlo(2), zhi(2)
           do i = zlo(1), zhi(1)
              Bz(i,k,1) = sigx1(i)*Bz(i,k,1) - sigx2(i)*(Ey(i+1,k,1)+Ey(i+1,k,2)+Ey(i+1,k,3) &
                   &                                    -Ey(i  ,k,1)-Ey(i  ,k,2)-Ey(i  ,k,3))
           end do
        end do

    else

        ! Cole-Karkkainen-Cowan push

        ! computes coefficients according to Cowan - PRST-AB 16, 041303 (2013)
        delta = max(dtsdx,dtsdz)
        rx = (dtsdx/delta)**2
        rz = (dtsdz/delta)**2
        betaxz = 0.125*rz
        betazx = 0.125*rx
        alphax = 1. - 2.*betaxz
        alphaz = 1. - 2.*betazx


        do    k = xlo(2), xhi(2)
           do i = xlo(1), xhi(1)
              Bx(i,k,2) = sigz1(k)*Bx(i,k,2) + sigz2(k)*( alphaz*(Ey(i  ,k+1,1)+Ey(i  ,k+1,2)+Ey(i  ,k+1,3)  &
                                                                 -Ey(i  ,k  ,1)-Ey(i  ,k  ,2)-Ey(i  ,k  ,3)) &
                                                        + betazx*(Ey(i+1,k+1,1)+Ey(i+1,k+1,2)+Ey(i+1,k+1,3)  &
                                                                 -Ey(i+1,k  ,1)-Ey(i+1,k  ,2)-Ey(i+1,k  ,3)  &
                                                                 +Ey(i-1,k+1,1)+Ey(i-1,k+1,2)+Ey(i-1,k+1,3)  &
                                                                 -Ey(i-1,k  ,1)-Ey(i-1,k  ,2)-Ey(i-1,k  ,3))) 
           end do
        end do

        do    k = ylo(2), yhi(2)
           do i = ylo(1), yhi(1)
              By(i,k,1) = sigz1(k)*By(i,k,1) - sigz2(k)*( alphaz*(Ex(i  ,k+1,1)+Ex(i  ,k+1,2)+Ex(i  ,k+1,3)  &
                                                                 -Ex(i  ,k  ,1)-Ex(i  ,k  ,2)-Ex(i  ,k  ,3)) &
                                                        + betazx*(Ex(i+1,k+1,1)+Ex(i+1,k+1,2)+Ex(i+1,k+1,3)  &
                                                                 -Ex(i+1,k  ,1)-Ex(i+1,k  ,2)-Ex(i+1,k  ,3)  &
                                                                 +Ex(i-1,k+1,1)+Ex(i-1,k+1,2)+Ex(i-1,k+1,3)  &
                                                                 -Ex(i-1,k  ,1)-Ex(i-1,k  ,2)-Ex(i-1,k  ,3)))

              By(i,k,2) = sigx1(i)*By(i,k,2) + sigx2(i)*( alphax*(Ez(i+1,k  ,1)+Ez(i+1,k  ,2)+Ez(i+1,k  ,3)  & 
                                                                 -Ez(i  ,k  ,1)-Ez(i  ,k  ,2)-Ez(i  ,k  ,3)) &
                                                        + betaxz*(Ez(i+1,k+1,1)+Ez(i+1,k+1,2)+Ez(i+1,k+1,3)  & 
                                                                 -Ez(i  ,k+1,1)-Ez(i  ,k+1,2)-Ez(i  ,k+1,3)  &
                                                                 +Ez(i+1,k-1,1)+Ez(i+1,k-1,2)+Ez(i+1,k-1,3)  & 
                                                                 -Ez(i  ,k-1,1)-Ez(i  ,k-1,2)-Ez(i  ,k-1,3)))
                   
           end do
        end do

        do    k = zlo(2), zhi(2)
           do i = zlo(1), zhi(1)
              Bz(i,k,1) = sigx1(i)*Bz(i,k,1) - sigx2(i)*( alphax*(Ey(i+1,k  ,1)+Ey(i+1,k  ,2)+Ey(i+1,k  ,3)  &
                                                                 -Ey(i  ,k  ,1)-Ey(i  ,k  ,2)-Ey(i  ,k  ,3)) &
                                                        + betaxz*(Ey(i+1,k+1,1)+Ey(i+1,k+1,2)+Ey(i+1,k+1,3)  &
                                                                 -Ey(i  ,k+1,1)-Ey(i  ,k+1,2)-Ey(i  ,k+1,3)  &
                                                                 +Ey(i+1,k-1,1)+Ey(i+1,k-1,2)+Ey(i+1,k-1,3)  &
                                                                 -Ey(i  ,k-1,1)-Ey(i  ,k-1,2)-Ey(i  ,k-1,3)))
           end do
        end do

    end if


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


  subroutine warpx_push_pml_f_3d (lo, hi, &
       &                          f ,  flo,  fhi, &
       &                          Ex, Exlo, Exhi, &
       &                          Ey, Eylo, Eyhi, &
       &                          Ez, Ezlo, Ezhi, &
       &                          sigx1, sigx1_lo, sigx1_hi, &
       &                          sigx2, sigx2_lo, sigx2_hi, &
       &                          sigy1, sigy1_lo, sigy1_hi, &
       &                          sigy2, sigy2_lo, sigy2_hi, &
       &                          sigz1, sigz1_lo, sigz1_hi, &
       &                          sigz2, sigz2_lo, sigz2_hi, c2inv) &
       bind(c,name='warpx_push_pml_f_3d')
    integer, intent(in) :: lo(3), hi(3), Exlo(3), Exhi(3), Eylo(3), Eyhi(3), Ezlo(3), Ezhi(3), &
         flo(3), fhi(3)
    integer, intent(in), value :: sigx1_lo, sigx1_hi, sigx2_lo, sigx2_hi
    integer, intent(in), value :: sigy1_lo, sigy1_hi, sigy2_lo, sigy2_hi
    integer, intent(in), value :: sigz1_lo, sigz1_hi, sigz2_lo, sigz2_hi
    real(amrex_real), intent(inout) :: f  ( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),3)
    real(amrex_real), intent(in   ) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),Exlo(3):Exhi(3),3)
    real(amrex_real), intent(in   ) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),Eylo(3):Eyhi(3),3)
    real(amrex_real), intent(in   ) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),Ezlo(3):Ezhi(3),3)
    real(amrex_real), intent(in) :: sigx1(sigx1_lo:sigx1_hi)
    real(amrex_real), intent(in) :: sigx2(sigx2_lo:sigx2_hi)
    real(amrex_real), intent(in) :: sigy1(sigy1_lo:sigy1_hi)
    real(amrex_real), intent(in) :: sigy2(sigy2_lo:sigy2_hi)
    real(amrex_real), intent(in) :: sigz1(sigz1_lo:sigz1_hi)
    real(amrex_real), intent(in) :: sigz2(sigz2_lo:sigz2_hi)
    real(amrex_real), intent(in) :: c2inv

    integer :: i, j, k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             f(i,j,k,1) = sigx1(i)*f(i,j,k,1) + sigx2(i)*c2inv* &
                  &                             ((Ex(i,j,k,1)-Ex(i-1,j,k,1)) &
                  &                            + (Ex(i,j,k,2)-Ex(i-1,j,k,2)) &
                  &                            + (Ex(i,j,k,3)-Ex(i-1,j,k,3)))
          end do
          do i = lo(1), hi(1)
             f(i,j,k,2) = sigy1(j)*f(i,j,k,2) + sigy2(j)*c2inv* &
                  &                             ((Ey(i,j,k,1)-Ey(i,j-1,k,1)) &
                  &                            + (Ey(i,j,k,2)-Ey(i,j-1,k,2)) &
                  &                            + (Ey(i,j,k,3)-Ey(i,j-1,k,3)))
          end do
          do i = lo(1), hi(1)
             f(i,j,k,3) = sigz1(k)*f(i,j,k,3) + sigz2(k)*c2inv* &
                  &                             ((Ez(i,j,k,1)-Ez(i,j,k-1,1)) &
                  &                            + (Ez(i,j,k,2)-Ez(i,j,k-1,2)) &
                  &                            + (Ez(i,j,k,3)-Ez(i,j,k-1,3)))
          end do
       end do
    end do
  end subroutine warpx_push_pml_f_3d

  subroutine warpx_push_pml_f_2d (lo, hi, &
       &                          f ,  flo,  fhi, &
       &                          Ex, Exlo, Exhi, &
       &                          Ey, Eylo, Eyhi, &
       &                          Ez, Ezlo, Ezhi, &
       &                          sigx1, sigx1_lo, sigx1_hi, &
       &                          sigx2, sigx2_lo, sigx2_hi, &
       &                          sigz1, sigz1_lo, sigz1_hi, &
       &                          sigz2, sigz2_lo, sigz2_hi, c2inv) &
       bind(c,name='warpx_push_pml_f_2d')
    integer, intent(in) :: lo(2), hi(2), Exlo(2), Exhi(2), Eylo(2), Eyhi(2), Ezlo(2), Ezhi(2), &
         flo(2), fhi(2)
    integer, intent(in), value :: sigx1_lo, sigx1_hi, sigx2_lo, sigx2_hi
    integer, intent(in), value :: sigz1_lo, sigz1_hi, sigz2_lo, sigz2_hi
    real(amrex_real), intent(inout) :: f  ( flo(1): fhi(1), flo(2): fhi(2),3)
    real(amrex_real), intent(in   ) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),3)
    real(amrex_real), intent(in   ) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),3)
    real(amrex_real), intent(in   ) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),3)
    real(amrex_real), intent(in) :: sigx1(sigx1_lo:sigx1_hi)
    real(amrex_real), intent(in) :: sigx2(sigx2_lo:sigx2_hi)
    real(amrex_real), intent(in) :: sigz1(sigz1_lo:sigz1_hi)
    real(amrex_real), intent(in) :: sigz2(sigz2_lo:sigz2_hi)
    real(amrex_real), intent(in) :: c2inv

    integer :: i, k

    do    k = lo(2), hi(2)
       do i = lo(1), hi(1)
          f(i,k,1) = sigx1(i)*f(i,k,1) + sigx2(i)*c2inv* &
               &                         ((Ex(i,k,1)-Ex(i-1,k,1)) &
               &                        + (Ex(i,k,2)-Ex(i-1,k,2)) &
               &                        + (Ex(i,k,3)-Ex(i-1,k,3)))
       end do
       do i = lo(1), hi(1)
          f(i,k,3) = sigz1(k)*f(i,k,3) + sigz2(k)*c2inv* &
               &                         ((Ez(i,k,1)-Ez(i,k-1,1)) &
               &                        + (Ez(i,k,2)-Ez(i,k-1,2)) &
               &                        + (Ez(i,k,3)-Ez(i,k-1,3)))
       end do
    end do
  end subroutine warpx_push_pml_f_2d


  subroutine warpx_push_pml_evec_f_3d (xlo, xhi, ylo, yhi, zlo, zhi, &
       &                             Ex, Exlo, Exhi, &
       &                             Ey, Eylo, Eyhi, &
       &                             Ez, Ezlo, Ezhi, &
       &                              f,  flo,  fhi, &
       &                             sigx1, sigx1_lo, sigx1_hi, &
       &                             sigx2, sigx2_lo, sigx2_hi, &
       &                             sigy1, sigy1_lo, sigy1_hi, &
       &                             sigy2, sigy2_lo, sigy2_hi, &
       &                             sigz1, sigz1_lo, sigz1_hi, &
       &                             sigz2, sigz2_lo, sigz2_hi, c2) &
       bind(c,name='warpx_push_pml_evec_f_3d')
    integer, intent(in) :: xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
         Exlo(3), Exhi(3), Eylo(3), Eyhi(3), Ezlo(3), Ezhi(3), flo(3), fhi(3)
    integer, intent(in), value :: sigx1_lo, sigx1_hi, sigx2_lo, sigx2_hi
    integer, intent(in), value :: sigy1_lo, sigy1_hi, sigy2_lo, sigy2_hi
    integer, intent(in), value :: sigz1_lo, sigz1_hi, sigz2_lo, sigz2_hi
    real(amrex_real), intent(inout) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),Exlo(3):Exhi(3),3)
    real(amrex_real), intent(inout) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),Eylo(3):Eyhi(3),3)
    real(amrex_real), intent(inout) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),Ezlo(3):Ezhi(3),3)
    real(amrex_real), intent(in   ) ::  f ( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),3)
    real(amrex_real), intent(in) :: sigx1(sigx1_lo:sigx1_hi)
    real(amrex_real), intent(in) :: sigx2(sigx2_lo:sigx2_hi)
    real(amrex_real), intent(in) :: sigy1(sigy1_lo:sigy1_hi)
    real(amrex_real), intent(in) :: sigy2(sigy2_lo:sigy2_hi)
    real(amrex_real), intent(in) :: sigz1(sigz1_lo:sigz1_hi)
    real(amrex_real), intent(in) :: sigz2(sigz2_lo:sigz2_hi)
    real(amrex_real), intent(in) :: c2

    integer :: i, j, k

    do       k = xlo(3), xhi(3)
       do    j = xlo(2), xhi(2)
          do i = xlo(1), xhi(1)
             Ex(i,j,k,3) = sigx1(i)*Ex(i,j,k,3) + sigx2(i)*c2*((f(i+1,j,k,1)-f(i,j,k,1)) &
                  &                                          + (f(i+1,j,k,2)-f(i,j,k,2)) &
                  &                                          + (f(i+1,j,k,3)-f(i,j,k,3)))
          end do
       end do
    end do

    do       k = ylo(3), yhi(3)
       do    j = ylo(2), yhi(2)
          do i = ylo(1), yhi(1)
             Ey(i,j,k,3) = sigy1(j)*Ey(i,j,k,3) + sigy2(j)*c2*((f(i,j+1,k,1)-f(i,j,k,1)) &
                  &                                          + (f(i,j+1,k,2)-f(i,j,k,2)) &
                  &                                          + (f(i,j+1,k,3)-f(i,j,k,3)))
          end do
       end do
    end do

    do       k = zlo(3), zhi(3)
       do    j = zlo(2), zhi(2)
          do i = zlo(1), zhi(1)
             Ez(i,j,k,3) = sigz1(k)*Ez(i,j,k,3) + sigz2(k)*c2*((f(i,j,k+1,1)-f(i,j,k,1)) &
                  &                                          + (f(i,j,k+1,2)-f(i,j,k,2)) &
                  &                                          + (f(i,j,k+1,3)-f(i,j,k,3)))
          end do
       end do
    end do

  end subroutine warpx_push_pml_evec_f_3d


  subroutine warpx_push_pml_evec_f_2d (xlo, xhi, ylo, yhi, zlo, zhi, &
       &                             Ex, Exlo, Exhi, &
       &                             Ey, Eylo, Eyhi, &
       &                             Ez, Ezlo, Ezhi, &
       &                              f,  flo,  fhi, &
       &                             sigx1, sigx1_lo, sigx1_hi, &
       &                             sigx2, sigx2_lo, sigx2_hi, &
       &                             sigz1, sigz1_lo, sigz1_hi, &
       &                             sigz2, sigz2_lo, sigz2_hi, c2) &
       bind(c,name='warpx_push_pml_evec_f_2d')
    integer, intent(in) :: xlo(2), xhi(2), ylo(2), yhi(2), zlo(2), zhi(2), &
         Exlo(2), Exhi(2), Eylo(2), Eyhi(2), Ezlo(2), Ezhi(2), flo(2), fhi(2)
    integer, intent(in), value :: sigx1_lo, sigx1_hi, sigx2_lo, sigx2_hi
    integer, intent(in), value :: sigz1_lo, sigz1_hi, sigz2_lo, sigz2_hi
    real(amrex_real), intent(inout) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),3)
    real(amrex_real), intent(inout) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),3)
    real(amrex_real), intent(inout) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),3)
    real(amrex_real), intent(in   ) ::  f ( flo(1): fhi(1), flo(2): fhi(2),3)
    real(amrex_real), intent(in) :: sigx1(sigx1_lo:sigx1_hi)
    real(amrex_real), intent(in) :: sigx2(sigx2_lo:sigx2_hi)
    real(amrex_real), intent(in) :: sigz1(sigz1_lo:sigz1_hi)
    real(amrex_real), intent(in) :: sigz2(sigz2_lo:sigz2_hi)
    real(amrex_real), intent(in) :: c2

    integer :: i, k

    do    k = xlo(2), xhi(2)
       do i = xlo(1), xhi(1)
          Ex(i,k,3) = sigx1(i)*Ex(i,k,3) + sigx2(i)*c2*((f(i+1,k,1)-f(i,k,1)) &
               &                                      + (f(i+1,k,2)-f(i,k,2)) &
               &                                      + (f(i+1,k,3)-f(i,k,3)))
       end do
    end do

    do    k = zlo(2), zhi(2)
       do i = zlo(1), zhi(1)
          Ez(i,k,3) = sigz1(k)*Ez(i,k,3) + sigz2(k)*c2*((f(i,k+1,1)-f(i,k,1)) &
               &                                      + (f(i,k+1,2)-f(i,k,2)) &
               &                                      + (f(i,k+1,3)-f(i,k,3)))
       end do
    end do

  end subroutine warpx_push_pml_evec_f_2d

end module warpx_pml_module
