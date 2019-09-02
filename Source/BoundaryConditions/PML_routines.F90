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
       &                             dtsdx, dtsdy, dtsdz, solver_type) &
       bind(c,name='warpx_push_pml_bvec_3d')
    use amrex_constants_module, only : one, two, four, eighth
    integer, intent(in) :: xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
         Exlo(3), Exhi(3), Eylo(3), Eyhi(3), Ezlo(3), Ezhi(3), &
         Bxlo(3), Bxhi(3), Bylo(3), Byhi(3), Bzlo(3), Bzhi(3), solver_type
    real(amrex_real), intent(in   ) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),Exlo(3):Exhi(3),3)
    real(amrex_real), intent(in   ) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),Eylo(3):Eyhi(3),3)
    real(amrex_real), intent(in   ) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),Ezlo(3):Ezhi(3),3)
    real(amrex_real), intent(inout) :: Bx (Bxlo(1):Bxhi(1),Bxlo(2):Bxhi(2),Bxlo(3):Bxhi(3),2)
    real(amrex_real), intent(inout) :: By (Bylo(1):Byhi(1),Bylo(2):Byhi(2),Bylo(3):Byhi(3),2)
    real(amrex_real), intent(inout) :: Bz (Bzlo(1):Bzhi(1),Bzlo(2):Bzhi(2),Bzlo(3):Bzhi(3),2)
    real(amrex_real), intent(in) :: dtsdx, dtsdy, dtsdz

    real(amrex_real), parameter :: sixteenth = 1.d0/16.d0

    integer :: i, j, k

    real(amrex_real) :: delta, rx, ry, rz, betaxz, betaxy, betayx, betayz, betazx, betazy
    real(amrex_real) :: beta, alphax, alphay, alphaz, gammax, gammay, gammaz

    ! solver_type: 0=Yee; 1=CKC

    if (solver_type==0) then

        ! Yee push

        do       k = xlo(3), xhi(3)
           do    j = xlo(2), xhi(2)
              do i = xlo(1), xhi(1)
                 Bx(i,j,k,1) = Bx(i,j,k,1) - dtsdy*(Ez(i,j+1,k  ,1)+Ez(i,j+1,k  ,2)+Ez(i,j+1,k  ,3) &
                      &                            -Ez(i,j  ,k  ,1)-Ez(i,j  ,k  ,2)-Ez(i,j  ,k  ,3))
                 Bx(i,j,k,2) = Bx(i,j,k,2) + dtsdz*(Ey(i,j  ,k+1,1)+Ey(i,j  ,k+1,2)+Ey(i,j  ,k+1,3) &
                      &                            -Ey(i,j  ,k  ,1)-Ey(i,j  ,k  ,2)-Ey(i,j  ,k  ,3))
              end do
           end do
        end do

        do       k = ylo(3), yhi(3)
           do    j = ylo(2), yhi(2)
              do i = ylo(1), yhi(1)
                 By(i,j,k,1) = By(i,j,k,1) - dtsdz*(Ex(i  ,j,k+1,1)+Ex(i  ,j,k+1,2)+Ex(i  ,j,k+1,3) &
                      &                            -Ex(i  ,j,k  ,1)-Ex(i  ,j,k  ,2)-Ex(i  ,j,k  ,3))
                 By(i,j,k,2) = By(i,j,k,2) + dtsdx*(Ez(i+1,j,k  ,1)+Ez(i+1,j,k  ,2)+Ez(i+1,j,k  ,3) &
                      &                            -Ez(i  ,j,k  ,1)-Ez(i  ,j,k  ,2)-Ez(i  ,j,k  ,3))
              end do
           end do
        end do

        do       k = zlo(3), zhi(3)
           do    j = zlo(2), zhi(2)
              do i = zlo(1), zhi(1)
                 Bz(i,j,k,1) = Bz(i,j,k,1) - dtsdx*(Ey(i+1,j  ,k,1)+Ey(i+1,j  ,k,2)+Ey(i+1,j  ,k,3) &
                      &                            -Ey(i  ,j  ,k,1)-Ey(i  ,j  ,k,2)-Ey(i  ,j  ,k,3))
                 Bz(i,j,k,2) = Bz(i,j,k,2) + dtsdy*(Ex(i  ,j+1,k,1)+Ex(i  ,j+1,k,2)+Ex(i  ,j+1,k,3) &
                      &                            -Ex(i  ,j  ,k,1)-Ex(i  ,j  ,k,2)-Ex(i  ,j  ,k,3))
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
        beta = eighth*(one-rx*ry*rz/(ry*rz+rz*rx+rx*ry))
        betaxy = ry*beta
        betaxz = rz*beta
        betayx = rx*beta
        betayz = rz*beta
        betazx = rx*beta
        betazy = ry*beta
        gammax = ry*rz*(sixteenth-eighth*ry*rz/(ry*rz+rz*rx+rx*ry))
        gammay = rx*rz*(sixteenth-eighth*rx*rz/(ry*rz+rz*rx+rx*ry))
        gammaz = rx*ry*(sixteenth-eighth*rx*ry/(ry*rz+rz*rx+rx*ry))
        alphax = one - two*betaxy - two* betaxz - four*gammax
        alphay = one - two*betayx - two* betayz - four*gammay
        alphaz = one - two*betazx - two* betazy - four*gammaz

        betaxy = dtsdx*betaxy
        betaxz = dtsdx*betaxz
        betayx = dtsdy*betayx
        betayz = dtsdy*betayz
        betazx = dtsdz*betazx
        betazy = dtsdz*betazy
        alphax = dtsdx*alphax
        alphay = dtsdy*alphay
        alphaz = dtsdz*alphaz
        gammax = dtsdx*gammax
        gammay = dtsdy*gammay
        gammaz = dtsdz*gammaz

        do       k = xlo(3), xhi(3)
           do    j = xlo(2), xhi(2)
              do i = xlo(1), xhi(1)
                 Bx(i,j,k,1) = Bx(i,j,k,1) - (alphay*(Ez(i  ,j+1,k  ,1)+Ez(i  ,j+1,k  ,2)+Ez(i  ,j+1,k  ,3)  &
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


                 Bx(i,j,k,2) = Bx(i,j,k,2) + (alphaz*(Ey(i  ,j  ,k+1,1)+Ey(i  ,j  ,k+1,2)+Ey(i  ,j  ,k+1,3)  &
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
                 By(i,j,k,1) = By(i,j,k,1) - (alphaz*(Ex(i  ,j  ,k+1,1)+Ex(i  ,j  ,k+1,2)+Ex(i  ,j  ,k+1,3)  &
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


                 By(i,j,k,2) = By(i,j,k,2) + (alphax*(Ez(i+1,j  ,k  ,1)+Ez(i+1,j  ,k  ,2)+Ez(i+1,j  ,k  ,3)  &
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
                 Bz(i,j,k,1) = Bz(i,j,k,1) - (alphax*(Ey(i+1,j  ,k  ,1)+Ey(i+1,j  ,k  ,2)+Ey(i+1,j  ,k  ,3)  &
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

                 Bz(i,j,k,2) = Bz(i,j,k,2) + (alphay*(Ex(i  ,j+1,k  ,1)+Ex(i  ,j+1,k  ,2)+Ex(i  ,j+1,k  ,3)  &
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
       &                             dtsdx, dtsdy, dtsdz) &
       bind(c,name='warpx_push_pml_evec_3d')
    integer, intent(in) :: xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
         Exlo(3), Exhi(3), Eylo(3), Eyhi(3), Ezlo(3), Ezhi(3), &
         Bxlo(3), Bxhi(3), Bylo(3), Byhi(3), Bzlo(3), Bzhi(3)
    real(amrex_real), intent(in) :: dtsdx, dtsdy, dtsdz
    real(amrex_real), intent(inout) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),Exlo(3):Exhi(3),2)
    real(amrex_real), intent(inout) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),Eylo(3):Eyhi(3),2)
    real(amrex_real), intent(inout) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),Ezlo(3):Ezhi(3),2)
    real(amrex_real), intent(in   ) :: Bx (Bxlo(1):Bxhi(1),Bxlo(2):Bxhi(2),Bxlo(3):Bxhi(3),2)
    real(amrex_real), intent(in   ) :: By (Bylo(1):Byhi(1),Bylo(2):Byhi(2),Bylo(3):Byhi(3),2)
    real(amrex_real), intent(in   ) :: Bz (Bzlo(1):Bzhi(1),Bzlo(2):Bzhi(2),Bzlo(3):Bzhi(3),2)

    integer :: i, j, k

    do       k = xlo(3), xhi(3)
       do    j = xlo(2), xhi(2)
          do i = xlo(1), xhi(1)
             Ex(i,j,k,1) = Ex(i,j,k,1) + dtsdy*(Bz(i,j  ,k  ,1)+Bz(i,j  ,k  ,2) &
                  &                            -Bz(i,j-1,k  ,1)-Bz(i,j-1,k  ,2))
             Ex(i,j,k,2) = Ex(i,j,k,2) - dtsdz*(By(i,j  ,k  ,1)+By(i,j  ,k  ,2) &
                  &                            -By(i,j  ,k-1,1)-By(i,j  ,k-1,2))
          end do
       end do
    end do

    do       k = ylo(3), yhi(3)
       do    j = ylo(2), yhi(2)
          do i = ylo(1), yhi(1)
             Ey(i,j,k,1) = Ey(i,j,k,1) + dtsdz*(Bx(i  ,j,k  ,1)+Bx(i  ,j,k  ,2) &
                  &                            -Bx(i  ,j,k-1,1)-Bx(i  ,j,k-1,2))
             Ey(i,j,k,2) = Ey(i,j,k,2) - dtsdx*(Bz(i  ,j,k  ,1)+Bz(i  ,j,k  ,2) &
                  &                            -Bz(i-1,j,k  ,1)-Bz(i-1,j,k  ,2))
          end do
       end do
    end do

    do       k = zlo(3), zhi(3)
       do    j = zlo(2), zhi(2)
          do i = zlo(1), zhi(1)
             Ez(i,j,k,1) = Ez(i,j,k,1) + dtsdx*(By(i  ,j  ,k,1)+By(i  ,j  ,k,2) &
                  &                            -By(i-1,j  ,k,1)-By(i-1,j  ,k,2))
             Ez(i,j,k,2) = Ez(i,j,k,2) - dtsdy*(Bx(i  ,j  ,k,1)+Bx(i  ,j  ,k,2) &
                  &                            -Bx(i  ,j-1,k,1)-Bx(i  ,j-1,k,2))
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
       &                             dtsdx, dtsdy, dtsdz, solver_type) &
       bind(c,name='warpx_push_pml_bvec_2d')
    use amrex_constants_module, only : one, two, eighth
    integer, intent(in) :: xlo(2), xhi(2), ylo(2), yhi(2), zlo(2), zhi(2), &
         Exlo(2), Exhi(2), Eylo(2), Eyhi(2), Ezlo(2), Ezhi(2), &
         Bxlo(2), Bxhi(2), Bylo(2), Byhi(2), Bzlo(2), Bzhi(2), solver_type
    real(amrex_real), intent(in   ) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),3)
    real(amrex_real), intent(in   ) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),3)
    real(amrex_real), intent(in   ) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),3)
    real(amrex_real), intent(inout) :: Bx (Bxlo(1):Bxhi(1),Bxlo(2):Bxhi(2),2)
    real(amrex_real), intent(inout) :: By (Bylo(1):Byhi(1),Bylo(2):Byhi(2),2)
    real(amrex_real), intent(inout) :: Bz (Bzlo(1):Bzhi(1),Bzlo(2):Bzhi(2),2)
    real(amrex_real), intent(in) :: dtsdx, dtsdy, dtsdz

    integer :: i, k

    real(amrex_real) :: delta, rx, rz, betaxz, betazx, alphax, alphaz


    ! solver_type: 0=Yee; 1=CKC

    if (solver_type==0) then

    ! Yee push

        do    k = xlo(2), xhi(2)
           do i = xlo(1), xhi(1)
              Bx(i,k,2) = Bx(i,k,2) + dtsdz*(Ey(i,k+1,1)+Ey(i,k+1,2)+Ey(i,k+1,3) &
                   &                        -Ey(i,k  ,1)-Ey(i,k  ,2)-Ey(i,k  ,3))
           end do
        end do

        do    k = ylo(2), yhi(2)
           do i = ylo(1), yhi(1)
              By(i,k,1) = By(i,k,1) - dtsdz*(Ex(i  ,k+1,1)+Ex(i  ,k+1,2)+Ex(i  ,k+1,3) &
                   &                        -Ex(i  ,k  ,1)-Ex(i  ,k  ,2)-Ex(i  ,k  ,3))
              By(i,k,2) = By(i,k,2) + dtsdx*(Ez(i+1,k  ,1)+Ez(i+1,k  ,2)+Ez(i+1,k  ,3) &
                   &                        -Ez(i  ,k  ,1)-Ez(i  ,k  ,2)-Ez(i  ,k  ,3))
           end do
        end do

        do    k = zlo(2), zhi(2)
           do i = zlo(1), zhi(1)
              Bz(i,k,1) = Bz(i,k,1) - dtsdx*(Ey(i+1,k,1)+Ey(i+1,k,2)+Ey(i+1,k,3) &
                   &                        -Ey(i  ,k,1)-Ey(i  ,k,2)-Ey(i  ,k,3))
           end do
        end do

    else

        ! Cole-Karkkainen-Cowan push

        ! computes coefficients according to Cowan - PRST-AB 16, 041303 (2013)
        delta = max(dtsdx,dtsdz)
        rx = (dtsdx/delta)**2
        rz = (dtsdz/delta)**2
        betaxz = eighth*rz
        betazx = eighth*rx
        alphax = one - two*betaxz
        alphaz = one - two*betazx

        betaxz = dtsdx*betaxz
        betazx = dtsdz*betazx
        alphax = dtsdx*alphax
        alphaz = dtsdz*alphaz

        do    k = xlo(2), xhi(2)
           do i = xlo(1), xhi(1)
              Bx(i,k,2) = Bx(i,k,2) + ( alphaz*(Ey(i  ,k+1,1)+Ey(i  ,k+1,2)+Ey(i  ,k+1,3)  &
                                               -Ey(i  ,k  ,1)-Ey(i  ,k  ,2)-Ey(i  ,k  ,3)) &
                                      + betazx*(Ey(i+1,k+1,1)+Ey(i+1,k+1,2)+Ey(i+1,k+1,3)  &
                                               -Ey(i+1,k  ,1)-Ey(i+1,k  ,2)-Ey(i+1,k  ,3)  &
                                               +Ey(i-1,k+1,1)+Ey(i-1,k+1,2)+Ey(i-1,k+1,3)  &
                                               -Ey(i-1,k  ,1)-Ey(i-1,k  ,2)-Ey(i-1,k  ,3)))
           end do
        end do

        do    k = ylo(2), yhi(2)
           do i = ylo(1), yhi(1)
              By(i,k,1) = By(i,k,1) - ( alphaz*(Ex(i  ,k+1,1)+Ex(i  ,k+1,2)+Ex(i  ,k+1,3)  &
                                               -Ex(i  ,k  ,1)-Ex(i  ,k  ,2)-Ex(i  ,k  ,3)) &
                                      + betazx*(Ex(i+1,k+1,1)+Ex(i+1,k+1,2)+Ex(i+1,k+1,3)  &
                                               -Ex(i+1,k  ,1)-Ex(i+1,k  ,2)-Ex(i+1,k  ,3)  &
                                               +Ex(i-1,k+1,1)+Ex(i-1,k+1,2)+Ex(i-1,k+1,3)  &
                                               -Ex(i-1,k  ,1)-Ex(i-1,k  ,2)-Ex(i-1,k  ,3)))

              By(i,k,2) = By(i,k,2) + ( alphax*(Ez(i+1,k  ,1)+Ez(i+1,k  ,2)+Ez(i+1,k  ,3)  &
                                               -Ez(i  ,k  ,1)-Ez(i  ,k  ,2)-Ez(i  ,k  ,3)) &
                                      + betaxz*(Ez(i+1,k+1,1)+Ez(i+1,k+1,2)+Ez(i+1,k+1,3)  &
                                               -Ez(i  ,k+1,1)-Ez(i  ,k+1,2)-Ez(i  ,k+1,3)  &
                                               +Ez(i+1,k-1,1)+Ez(i+1,k-1,2)+Ez(i+1,k-1,3)  &
                                               -Ez(i  ,k-1,1)-Ez(i  ,k-1,2)-Ez(i  ,k-1,3)))

           end do
        end do

        do    k = zlo(2), zhi(2)
           do i = zlo(1), zhi(1)
              Bz(i,k,1) = Bz(i,k,1) - ( alphax*(Ey(i+1,k  ,1)+Ey(i+1,k  ,2)+Ey(i+1,k  ,3)  &
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
       &                             dtsdx, dtsdy, dtsdz) &
       bind(c,name='warpx_push_pml_evec_2d')
    integer, intent(in) :: xlo(2), xhi(2), ylo(2), yhi(2), zlo(2), zhi(2), &
         Exlo(2), Exhi(2), Eylo(2), Eyhi(2), Ezlo(2), Ezhi(2), &
         Bxlo(2), Bxhi(2), Bylo(2), Byhi(2), Bzlo(2), Bzhi(2)
    real(amrex_real), intent(in) :: dtsdx, dtsdy, dtsdz
    real(amrex_real), intent(inout) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),2)
    real(amrex_real), intent(inout) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),2)
    real(amrex_real), intent(inout) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),2)
    real(amrex_real), intent(in   ) :: Bx (Bxlo(1):Bxhi(1),Bxlo(2):Bxhi(2),2)
    real(amrex_real), intent(in   ) :: By (Bylo(1):Byhi(1),Bylo(2):Byhi(2),2)
    real(amrex_real), intent(in   ) :: Bz (Bzlo(1):Bzhi(1),Bzlo(2):Bzhi(2),2)

    integer :: i, k

    do    k = xlo(2), xhi(2)
       do i = xlo(1), xhi(1)
          Ex(i,k,2) = Ex(i,k,2) - dtsdz*(By(i,k  ,1)+By(i,k  ,2) &
               &                        -By(i,k-1,1)-By(i,k-1,2))
       end do
    end do

    do    k = ylo(2), yhi(2)
       do i = ylo(1), yhi(1)
          Ey(i,k,1) = Ey(i,k,1) + dtsdz*(Bx(i  ,k  ,1)+Bx(i  ,k  ,2) &
               &                        -Bx(i  ,k-1,1)-Bx(i  ,k-1,2))
          Ey(i,k,2) = Ey(i,k,2) - dtsdx*(Bz(i  ,k  ,1)+Bz(i  ,k  ,2) &
               &                        -Bz(i-1,k  ,1)-Bz(i-1,k  ,2))
       end do
    end do

    do    k = zlo(2), zhi(2)
       do i = zlo(1), zhi(1)
          Ez(i,k,1) = Ez(i,k,1) + dtsdx*(By(i  ,k,1)+By(i  ,k,2) &
               &                        -By(i-1,k,1)-By(i-1,k,2))
       end do
    end do

  end subroutine warpx_push_pml_evec_2d


  subroutine warpx_push_pml_f_3d (lo, hi, &
       &                          f ,  flo,  fhi, &
       &                          Ex, Exlo, Exhi, &
       &                          Ey, Eylo, Eyhi, &
       &                          Ez, Ezlo, Ezhi, &
       &                          dtdx, dtdy, dtdz) &
       bind(c,name='warpx_push_pml_f_3d')
    integer, intent(in) :: lo(3), hi(3), Exlo(3), Exhi(3), Eylo(3), Eyhi(3), Ezlo(3), Ezhi(3), &
         flo(3), fhi(3)
    real(amrex_real), intent(inout) :: f  ( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),3)
    real(amrex_real), intent(in   ) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),Exlo(3):Exhi(3),3)
    real(amrex_real), intent(in   ) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),Eylo(3):Eyhi(3),3)
    real(amrex_real), intent(in   ) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),Ezlo(3):Ezhi(3),3)
    real(amrex_real), intent(in) :: dtdx, dtdy, dtdz

    integer :: i, j, k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             f(i,j,k,1) = f(i,j,k,1) + dtdx*((Ex(i,j,k,1)-Ex(i-1,j,k,1)) &
                  &                        + (Ex(i,j,k,2)-Ex(i-1,j,k,2)) &
                  &                        + (Ex(i,j,k,3)-Ex(i-1,j,k,3)))
             f(i,j,k,2) = f(i,j,k,2) + dtdy*((Ey(i,j,k,1)-Ey(i,j-1,k,1)) &
                  &                        + (Ey(i,j,k,2)-Ey(i,j-1,k,2)) &
                  &                        + (Ey(i,j,k,3)-Ey(i,j-1,k,3)))
             f(i,j,k,3) = f(i,j,k,3) + dtdz*((Ez(i,j,k,1)-Ez(i,j,k-1,1)) &
                  &                        + (Ez(i,j,k,2)-Ez(i,j,k-1,2)) &
                  &                        + (Ez(i,j,k,3)-Ez(i,j,k-1,3)))
          end do
       end do
    end do
  end subroutine warpx_push_pml_f_3d

  subroutine warpx_push_pml_f_2d (lo, hi, &
       &                          f ,  flo,  fhi, &
       &                          Ex, Exlo, Exhi, &
       &                          Ey, Eylo, Eyhi, &
       &                          Ez, Ezlo, Ezhi, &
       &                          dtdx, dtdy, dtdz) &
       bind(c,name='warpx_push_pml_f_2d')
    integer, intent(in) :: lo(2), hi(2), Exlo(2), Exhi(2), Eylo(2), Eyhi(2), Ezlo(2), Ezhi(2), &
         flo(2), fhi(2)
    real(amrex_real), intent(inout) :: f  ( flo(1): fhi(1), flo(2): fhi(2),3)
    real(amrex_real), intent(in   ) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),3)
    real(amrex_real), intent(in   ) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),3)
    real(amrex_real), intent(in   ) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),3)
    real(amrex_real), intent(in) :: dtdx, dtdy, dtdz

    integer :: i, k

    do    k = lo(2), hi(2)
       do i = lo(1), hi(1)
          f(i,k,1) = f(i,k,1) + dtdx*((Ex(i,k,1)-Ex(i-1,k,1)) &
               &                    + (Ex(i,k,2)-Ex(i-1,k,2)) &
               &                    + (Ex(i,k,3)-Ex(i-1,k,3)))
          f(i,k,3) = f(i,k,3) + dtdz*((Ez(i,k,1)-Ez(i,k-1,1)) &
               &                    + (Ez(i,k,2)-Ez(i,k-1,2)) &
               &                    + (Ez(i,k,3)-Ez(i,k-1,3)))
       end do
    end do
  end subroutine warpx_push_pml_f_2d

  subroutine warpx_push_pml_evec_f_3d (xlo, xhi, ylo, yhi, zlo, zhi, &
       &                             Ex, Exlo, Exhi, &
       &                             Ey, Eylo, Eyhi, &
       &                             Ez, Ezlo, Ezhi, &
       &                              f,  flo,  fhi, &
       &                             dtsdx, dtsdy, dtsdz, solver_type) &
       bind(c,name='warpx_push_pml_evec_f_3d')
    use amrex_constants_module, only : one, two, four, eighth
    integer, intent(in) :: xlo(3), xhi(3), ylo(3), yhi(3), zlo(3), zhi(3), &
         Exlo(3), Exhi(3), Eylo(3), Eyhi(3), Ezlo(3), Ezhi(3), flo(3), fhi(3), &
         solver_type
    real(amrex_real), intent(inout) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),Exlo(3):Exhi(3),3)
    real(amrex_real), intent(inout) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),Eylo(3):Eyhi(3),3)
    real(amrex_real), intent(inout) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),Ezlo(3):Ezhi(3),3)
    real(amrex_real), intent(in   ) ::  f ( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),3)
    real(amrex_real), intent(in) :: dtsdx, dtsdy, dtsdz

    real(amrex_real), parameter :: sixteenth = 1.d0/16.d0

    integer :: i, j, k

    real(amrex_real) :: delta, rx, ry, rz, betaxz, betaxy, betayx, betayz, betazx, betazy
    real(amrex_real) :: beta, alphax, alphay, alphaz, gammax, gammay, gammaz

    ! solver_type: 0=Yee; 1=CKC

    if (solver_type==0) then

        do       k = xlo(3), xhi(3)
           do    j = xlo(2), xhi(2)
              do i = xlo(1), xhi(1)
                 Ex(i,j,k,3) = Ex(i,j,k,3) + dtsdx*((f(i+1,j,k,1)-f(i,j,k,1)) &
                   &                                + (f(i+1,j,k,2)-f(i,j,k,2)) &
                   &                                + (f(i+1,j,k,3)-f(i,j,k,3)))
              end do
           end do
        end do

        do       k = ylo(3), yhi(3)
           do    j = ylo(2), yhi(2)
              do i = ylo(1), yhi(1)
                 Ey(i,j,k,3) = Ey(i,j,k,3) + dtsdx*((f(i,j+1,k,1)-f(i,j,k,1)) &
                   &                                + (f(i,j+1,k,2)-f(i,j,k,2)) &
                   &                                + (f(i,j+1,k,3)-f(i,j,k,3)))
              end do
           end do
        end do

        do       k = zlo(3), zhi(3)
           do    j = zlo(2), zhi(2)
              do i = zlo(1), zhi(1)
                 Ez(i,j,k,3) = Ez(i,j,k,3) + dtsdx*((f(i,j,k+1,1)-f(i,j,k,1)) &
                   &                                + (f(i,j,k+1,2)-f(i,j,k,2)) &
                   &                                + (f(i,j,k+1,3)-f(i,j,k,3)))
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
        beta = eighth*(one-rx*ry*rz/(ry*rz+rz*rx+rx*ry))
        betaxy = ry*beta
        betaxz = rz*beta
        betayx = rx*beta
        betayz = rz*beta
        betazx = rx*beta
        betazy = ry*beta
        gammax = ry*rz*(sixteenth-eighth*ry*rz/(ry*rz+rz*rx+rx*ry))
        gammay = rx*rz*(sixteenth-eighth*rx*rz/(ry*rz+rz*rx+rx*ry))
        gammaz = rx*ry*(sixteenth-eighth*rx*ry/(ry*rz+rz*rx+rx*ry))
        alphax = one - two*betaxy - two* betaxz - four*gammax
        alphay = one - two*betayx - two* betayz - four*gammay
        alphaz = one - two*betazx - two* betazy - four*gammaz

        betaxy = dtsdx*betaxy
        betaxz = dtsdx*betaxz
        betayx = dtsdy*betayx
        betayz = dtsdy*betayz
        betazx = dtsdz*betazx
        betazy = dtsdz*betazy
        alphax = dtsdx*alphax
        alphay = dtsdy*alphay
        alphaz = dtsdz*alphaz
        gammax = dtsdx*gammax
        gammay = dtsdy*gammay
        gammaz = dtsdz*gammaz

        do       k = xlo(3), xhi(3)
           do    j = xlo(2), xhi(2)
              do i = xlo(1), xhi(1)
                 Ex(i,j,k,3) = Ex(i,j,k,3) &
                     +alphax*(f(i+1,j  ,k  ,1)+f(i+1,j  ,k  ,2)+f(i+1,j  ,k  ,3)  &
                             -f(i  ,j  ,k  ,1)-f(i  ,j  ,k  ,2)-f(i  ,j  ,k  ,3)) &
                     +betaxy*(f(i+1,j+1,k  ,1)+f(i+1,j+1,k  ,2)+f(i+1,j+1,k  ,3)  &
                             -f(i  ,j+1,k  ,1)-f(i  ,j+1,k  ,2)-f(i  ,j+1,k  ,3)  &
                             +f(i+1,j-1,k  ,1)+f(i+1,j-1,k  ,2)+f(i+1,j-1,k  ,3)  &
                             -f(i  ,j-1,k  ,1)-f(i  ,j-1,k  ,2)-f(i  ,j-1,k  ,3)) &
                     +betaxz*(f(i+1,j  ,k+1,1)+f(i+1,j  ,k+1,2)+f(i+1,j  ,k+1,3)  &
                             -f(i  ,j  ,k+1,1)-f(i  ,j  ,k+1,2)-f(i  ,j  ,k+1,3)  &
                             +f(i+1,j  ,k-1,1)+f(i+1,j  ,k-1,2)+f(i+1,j  ,k-1,3)  &
                             -f(i  ,j  ,k-1,1)-f(i  ,j  ,k-1,2)-f(i  ,j  ,k-1,3)) &
                     +gammax*(f(i+1,j+1,k+1,1)+f(i+1,j+1,k+1,2)+f(i+1,j+1,k+1,3)  &
                             -f(i  ,j+1,k+1,1)-f(i  ,j+1,k+1,2)-f(i  ,j+1,k+1,3)  &
                             +f(i+1,j-1,k+1,1)+f(i+1,j-1,k+1,2)+f(i+1,j-1,k+1,3)  &
                             -f(i  ,j-1,k+1,1)-f(i  ,j-1,k+1,2)-f(i  ,j-1,k+1,3)  &
                             +f(i+1,j+1,k-1,1)+f(i+1,j+1,k-1,2)+f(i+1,j+1,k-1,3)  &
                             -f(i  ,j+1,k-1,1)-f(i  ,j+1,k-1,2)-f(i  ,j+1,k-1,3)  &
                             +f(i+1,j-1,k-1,1)+f(i+1,j-1,k-1,2)+f(i+1,j-1,k-1,3)  &
                             -f(i  ,j-1,k-1,1)-f(i  ,j-1,k-1,2)-f(i  ,j-1,k-1,3))
              end do
           end do
        end do

        do       k = ylo(3), yhi(3)
           do    j = ylo(2), yhi(2)
              do i = ylo(1), yhi(1)
                 Ey(i,j,k,3) = Ey(i,j,k,3) &
                     +alphay*(f(i  ,j+1,k  ,1)+f(i  ,j+1,k  ,2)+f(i  ,j+1,k  ,3)  &
                             -f(i  ,j  ,k  ,1)-f(i  ,j  ,k  ,2)-f(i  ,j  ,k  ,3)) &
                     +betayx*(f(i+1,j+1,k  ,1)+f(i+1,j+1,k  ,2)+f(i+1,j+1,k  ,3)  &
                             -f(i+1,j  ,k  ,1)-f(i+1,j  ,k  ,2)-f(i+1,j  ,k  ,3)  &
                             +f(i-1,j+1,k  ,1)+f(i-1,j+1,k  ,2)+f(i-1,j+1,k  ,3)  &
                             -f(i-1,j  ,k  ,1)-f(i-1,j  ,k  ,2)-f(i-1,j  ,k  ,3)) &
                     +betayz*(f(i  ,j+1,k+1,1)+f(i  ,j+1,k+1,2)+f(i  ,j+1,k+1,3)  &
                             -f(i  ,j  ,k+1,1)-f(i  ,j  ,k+1,2)-f(i  ,j  ,k+1,3)  &
                             +f(i  ,j+1,k-1,1)+f(i  ,j+1,k-1,2)+f(i  ,j+1,k-1,3)  &
                             -f(i  ,j  ,k-1,1)-f(i  ,j  ,k-1,2)-f(i  ,j  ,k-1,3)) &
                     +gammay*(f(i+1,j+1,k+1,1)+f(i+1,j+1,k+1,2)+f(i+1,j+1,k+1,3)  &
                             -f(i+1,j  ,k+1,1)-f(i+1,j  ,k+1,2)-f(i+1,j  ,k+1,3)  &
                             +f(i-1,j+1,k+1,1)+f(i-1,j+1,k+1,2)+f(i-1,j+1,k+1,3)  &
                             -f(i-1,j  ,k+1,1)-f(i-1,j  ,k+1,2)-f(i-1,j  ,k+1,3)  &
                             +f(i+1,j+1,k-1,1)+f(i+1,j+1,k-1,2)+f(i+1,j+1,k-1,3)  &
                             -f(i+1,j  ,k-1,1)-f(i+1,j  ,k-1,2)-f(i+1,j  ,k-1,3)  &
                             +f(i-1,j+1,k-1,1)+f(i-1,j+1,k-1,2)+f(i-1,j+1,k-1,3)  &
                             -f(i-1,j  ,k-1,1)-f(i-1,j  ,k-1,2)-f(i-1,j  ,k-1,3))
              end do
           end do
        end do

        do       k = zlo(3), zhi(3)
           do    j = zlo(2), zhi(2)
              do i = zlo(1), zhi(1)
                 Ez(i,j,k,3) = Ez(i,j,k,3) &
                     +alphaz*(f(i  ,j  ,k+1,1)+f(i  ,j  ,k+1,2)+f(i  ,j  ,k+1,3)  &
                             -f(i  ,j  ,k  ,1)-f(i  ,j  ,k  ,2)-f(i  ,j  ,k  ,3)) &
                     +betazx*(f(i+1,j  ,k+1,1)+f(i+1,j  ,k+1,2)+f(i+1,j  ,k+1,3)  &
                             -f(i+1,j  ,k  ,1)-f(i+1,j  ,k  ,2)-f(i+1,j  ,k  ,3)  &
                             +f(i-1,j  ,k+1,1)+f(i-1,j  ,k+1,2)+f(i-1,j  ,k+1,3)  &
                             -f(i-1,j  ,k  ,1)-f(i-1,j  ,k  ,2)-f(i-1,j  ,k  ,3)) &
                     +betazy*(f(i  ,j+1,k+1,1)+f(i  ,j+1,k+1,2)+f(i  ,j+1,k+1,3)  &
                             -f(i  ,j+1,k  ,1)-f(i  ,j+1,k  ,2)-f(i  ,j+1,k  ,3)  &
                             +f(i  ,j-1,k+1,1)+f(i  ,j-1,k+1,2)+f(i  ,j-1,k+1,3)  &
                             -f(i  ,j-1,k  ,1)-f(i  ,j-1,k  ,2)-f(i  ,j-1,k  ,3)) &
                     +gammaz*(f(i+1,j+1,k+1,1)+f(i+1,j+1,k+1,2)+f(i+1,j+1,k+1,3)  &
                             -f(i+1,j+1,k  ,1)-f(i+1,j+1,k  ,2)-f(i+1,j+1,k  ,3)  &
                             +f(i-1,j+1,k+1,1)+f(i-1,j+1,k+1,2)+f(i-1,j+1,k+1,3)  &
                             -f(i-1,j+1,k  ,1)-f(i-1,j+1,k  ,2)-f(i-1,j+1,k  ,3)  &
                             +f(i+1,j-1,k+1,1)+f(i+1,j-1,k+1,2)+f(i+1,j-1,k+1,3)  &
                             -f(i+1,j-1,k  ,1)-f(i+1,j-1,k  ,2)-f(i+1,j-1,k  ,3)  &
                             +f(i-1,j-1,k+1,1)+f(i-1,j-1,k+1,2)+f(i-1,j-1,k+1,3)  &
                             -f(i-1,j-1,k  ,1)-f(i-1,j-1,k  ,2)-f(i-1,j-1,k  ,3))
              end do
           end do
        end do

    endif

  end subroutine warpx_push_pml_evec_f_3d


  subroutine warpx_push_pml_evec_f_2d (xlo, xhi, ylo, yhi, zlo, zhi, &
       &                             Ex, Exlo, Exhi, &
       &                             Ey, Eylo, Eyhi, &
       &                             Ez, Ezlo, Ezhi, &
       &                              f,  flo,  fhi, &
       &                             dtsdx, dtsdy, dtsdz, solver_type) &
       bind(c,name='warpx_push_pml_evec_f_2d')
    use amrex_constants_module, only : one, two, four, eighth
    integer, intent(in) :: xlo(2), xhi(2), ylo(2), yhi(2), zlo(2), zhi(2), &
         Exlo(2), Exhi(2), Eylo(2), Eyhi(2), Ezlo(2), Ezhi(2), flo(2), fhi(2), &
         solver_type
    real(amrex_real), intent(inout) :: Ex (Exlo(1):Exhi(1),Exlo(2):Exhi(2),3)
    real(amrex_real), intent(inout) :: Ey (Eylo(1):Eyhi(1),Eylo(2):Eyhi(2),3)
    real(amrex_real), intent(inout) :: Ez (Ezlo(1):Ezhi(1),Ezlo(2):Ezhi(2),3)
    real(amrex_real), intent(in   ) ::  f ( flo(1): fhi(1), flo(2): fhi(2),3)
    real(amrex_real), intent(in) :: dtsdx, dtsdy, dtsdz

    integer :: i, k

    real(amrex_real) :: delta, rx, rz, betaxz, betazx, alphax, alphaz

    ! solver_type: 0=Yee; 1=CKC

    if (solver_type==0) then

        do    k = xlo(2), xhi(2)
           do i = xlo(1), xhi(1)
              Ex(i,k,3) = Ex(i,k,3) +   dtsdx*((f(i+1,k,1)-f(i,k,1)) &
                   &                         + (f(i+1,k,2)-f(i,k,2)) &
                   &                         + (f(i+1,k,3)-f(i,k,3)))
           end do
        end do

        do    k = zlo(2), zhi(2)
           do i = zlo(1), zhi(1)
              Ez(i,k,3) = Ez(i,k,3) +   dtsdz*((f(i,k+1,1)-f(i,k,1)) &
                   &                         + (f(i,k+1,2)-f(i,k,2)) &
                   &                         + (f(i,k+1,3)-f(i,k,3)))
           end do
        end do

    else

        ! Cole-Karkkainen-Cowan push

        ! computes coefficients according to Cowan - PRST-AB 16, 041303 (2013)
        delta = max(dtsdx,dtsdz)
        rx = (dtsdx/delta)**2
        rz = (dtsdz/delta)**2
        betaxz = eighth*rz
 	    betazx = eighth*rx
        alphax = one - two*betaxz
        alphaz = one - two*betazx

        betaxz = dtsdx*betaxz
        betazx = dtsdz*betazx
        alphax = dtsdx*alphax
        alphaz = dtsdz*alphaz

        do    k = xlo(2), xhi(2)
           do i = xlo(1), xhi(1)
              Ex(i,k,3) = Ex(i,k,3) &
                      + alphax*(f(i+1,k  ,1)+f(i+1,k  ,2)+f(i+1,k  ,3)  &
                               -f(i  ,k  ,1)-f(i  ,k  ,2)-f(i  ,k  ,3)) &
                      + betaxz*(f(i+1,k+1,1)+f(i+1,k+1,2)+f(i+1,k+1,3)  &
                               -f(i  ,k+1,1)-f(i  ,k+1,2)-f(i  ,k+1,3)  &
                               +f(i+1,k-1,1)+f(i+1,k-1,2)+f(i+1,k-1,3)  &
                               -f(i  ,k-1,1)-f(i  ,k-1,2)-f(i  ,k-1,3))
           end do
        end do

        do    k = zlo(2), zhi(2)
           do i = zlo(1), zhi(1)
              Ez(i,k,3) = Ez(i,k,3) &
                      + alphaz*(f(i  ,k+1,1)+f(i  ,k+1,2)+f(i  ,k+1,3)  &
                               -f(i  ,k  ,1)-f(i  ,k  ,2)-f(i  ,k  ,3)) &
                      + betazx*(f(i+1,k+1,1)+f(i+1,k+1,2)+f(i+1,k+1,3)  &
                               -f(i+1,k  ,1)-f(i+1,k  ,2)-f(i+1,k  ,3)  &
                               +f(i-1,k+1,1)+f(i-1,k+1,2)+f(i-1,k+1,3)  &
                               -f(i-1,k  ,1)-f(i-1,k  ,2)-f(i-1,k  ,3))
           end do
        end do

    endif

  end subroutine warpx_push_pml_evec_f_2d


  subroutine warpx_damp_pml_2d (texlo, texhi, teylo, teyhi, tezlo, tezhi, &
       &                        tbxlo, tbxhi, tbylo, tbyhi, tbzlo, tbzhi, &
       &                        ex, exlo, exhi, ey, eylo, eyhi, ez, ezlo, ezhi, &
       &                        bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, &
       &                        sigex, sexlo, sexhi, sigez, sezlo, sezhi, &
       &                        sigcx, scxlo, scxhi, sigcz, sczlo, sczhi) &
       bind(c,name='warpx_damp_pml_2d')
    integer, dimension(2), intent(in) :: texlo, texhi, teylo, teyhi, tezlo, tezhi, &
         tbxlo, tbxhi, tbylo, tbyhi, tbzlo, tbzhi, &
         exlo, exhi, eylo, eyhi, ezlo, ezhi, bxlo, bxhi, bylo, byhi, bzlo, bzhi
    integer, intent(in), value :: sexlo, sexhi, sezlo, sezhi, &
         &                        scxlo, scxhi, sczlo, sczhi
    real(amrex_real), intent(inout) :: ex(exlo(1):exhi(1),exlo(2):exhi(2),3)
    real(amrex_real), intent(inout) :: ey(eylo(1):eyhi(1),eylo(2):eyhi(2),3)
    real(amrex_real), intent(inout) :: ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2),3)
    real(amrex_real), intent(inout) :: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),2)
    real(amrex_real), intent(inout) :: by(bylo(1):byhi(1),bylo(2):byhi(2),2)
    real(amrex_real), intent(inout) :: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),2)
    real(amrex_real), intent(in) :: sigex(sexlo:sexhi)
    real(amrex_real), intent(in) :: sigez(sezlo:sezhi)
    real(amrex_real), intent(in) :: sigcx(scxlo:scxhi)
    real(amrex_real), intent(in) :: sigcz(sczlo:sczhi)

    integer :: i,k

    do    k = texlo(2), texhi(2)
       do i = texlo(1), texhi(1)
          ex(i,k,2) = ex(i,k,2) * sigez(k)
          ex(i,k,3) = ex(i,k,3) * sigcx(i)
       end do
    end do

    do    k = teylo(2), teyhi(2)
       do i = teylo(1), teyhi(1)
          ey(i,k,1) = ey(i,k,1) * sigez(k)
          ey(i,k,2) = ey(i,k,2) * sigex(i)
       end do
    end do

    do    k = tezlo(2), tezhi(2)
       do i = tezlo(1), tezhi(1)
          ez(i,k,1) = ez(i,k,1) * sigex(i)
          ez(i,k,3) = ez(i,k,3) * sigcz(k)
       end do
    end do

    do    k = tbxlo(2), tbxhi(2)
       do i = tbxlo(1), tbxhi(1)
          bx(i,k,2) = bx(i,k,2) * sigcz(k)
       end do
    end do

    do    k = tbylo(2), tbyhi(2)
       do i = tbylo(1), tbyhi(1)
          by(i,k,1) = by(i,k,1) * sigcz(k)
          by(i,k,2) = by(i,k,2) * sigcx(i)
       end do
    end do

    do    k = tbzlo(2), tbzhi(2)
       do i = tbzlo(1), tbzhi(1)
          bz(i,k,1) = bz(i,k,1) * sigcx(i)
       end do
    end do

  end subroutine warpx_damp_pml_2d


  subroutine warpx_damp_pml_3d (texlo, texhi, teylo, teyhi, tezlo, tezhi, &
       &                        tbxlo, tbxhi, tbylo, tbyhi, tbzlo, tbzhi, &
       &                        ex, exlo, exhi, ey, eylo, eyhi, ez, ezlo, ezhi, &
       &                        bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, &
       &                        sigex, sexlo, sexhi, sigey, seylo, seyhi, sigez, sezlo, sezhi, &
       &                        sigcx, scxlo, scxhi, sigcy, scylo, scyhi, sigcz, sczlo, sczhi) &
       bind(c,name='warpx_damp_pml_3d')
    integer, dimension(3), intent(in) :: texlo, texhi, teylo, teyhi, tezlo, tezhi, &
         tbxlo, tbxhi, tbylo, tbyhi, tbzlo, tbzhi, &
         exlo, exhi, eylo, eyhi, ezlo, ezhi, bxlo, bxhi, bylo, byhi, bzlo, bzhi
    integer, intent(in), value :: sexlo, sexhi, seylo, seyhi, sezlo, sezhi, &
         &                        scxlo, scxhi, scylo, scyhi, sczlo, sczhi
    real(amrex_real), intent(inout) :: ex(exlo(1):exhi(1),exlo(2):exhi(2),exlo(3):exhi(3),3)
    real(amrex_real), intent(inout) :: ey(eylo(1):eyhi(1),eylo(2):eyhi(2),eylo(3):eyhi(3),3)
    real(amrex_real), intent(inout) :: ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2),ezlo(3):ezhi(3),3)
    real(amrex_real), intent(inout) :: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3),2)
    real(amrex_real), intent(inout) :: by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3),2)
    real(amrex_real), intent(inout) :: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3),2)
    real(amrex_real), intent(in) :: sigex(sexlo:sexhi)
    real(amrex_real), intent(in) :: sigey(seylo:seyhi)
    real(amrex_real), intent(in) :: sigez(sezlo:sezhi)
    real(amrex_real), intent(in) :: sigcx(scxlo:scxhi)
    real(amrex_real), intent(in) :: sigcy(scylo:scyhi)
    real(amrex_real), intent(in) :: sigcz(sczlo:sczhi)

    integer :: i,j,k

    do       k = texlo(3), texhi(3)
       do    j = texlo(2), texhi(2)
          do i = texlo(1), texhi(1)
             ex(i,j,k,1) = ex(i,j,k,1) * sigey(j)
             ex(i,j,k,2) = ex(i,j,k,2) * sigez(k)
             ex(i,j,k,3) = ex(i,j,k,3) * sigcx(i)
          end do
       end do
    end do

    do       k = teylo(3), teyhi(3)
       do    j = teylo(2), teyhi(2)
          do i = teylo(1), teyhi(1)
             ey(i,j,k,1) = ey(i,j,k,1) * sigez(k)
             ey(i,j,k,2) = ey(i,j,k,2) * sigex(i)
             ey(i,j,k,3) = ey(i,j,k,3) * sigcy(j)
          end do
       end do
    end do

    do       k = tezlo(3), tezhi(3)
       do    j = tezlo(2), tezhi(2)
          do i = tezlo(1), tezhi(1)
             ez(i,j,k,1) = ez(i,j,k,1) * sigex(i)
             ez(i,j,k,2) = ez(i,j,k,2) * sigey(j)
             ez(i,j,k,3) = ez(i,j,k,3) * sigcz(k)
          end do
       end do
    end do

    do       k = tbxlo(3), tbxhi(3)
       do    j = tbxlo(2), tbxhi(2)
          do i = tbxlo(1), tbxhi(1)
             bx(i,j,k,1) = bx(i,j,k,1) * sigcy(j)
             bx(i,j,k,2) = bx(i,j,k,2) * sigcz(k)
          end do
       end do
    end do

    do       k = tbylo(3), tbyhi(3)
       do    j = tbylo(2), tbyhi(2)
          do i = tbylo(1), tbyhi(1)
             by(i,j,k,1) = by(i,j,k,1) * sigcz(k)
             by(i,j,k,2) = by(i,j,k,2) * sigcx(i)
          end do
       end do
    end do

    do       k = tbzlo(3), tbzhi(3)
       do    j = tbzlo(2), tbzhi(2)
          do i = tbzlo(1), tbzhi(1)
             bz(i,j,k,1) = bz(i,j,k,1) * sigcx(i)
             bz(i,j,k,2) = bz(i,j,k,2) * sigcy(j)
          end do
       end do
    end do

  end subroutine warpx_damp_pml_3d

  subroutine warpx_damp_pml_f_2d (tndlo, tndhi, f, flo, fhi,&
       &                          sigex, sexlo, sexhi, sigez, sezlo, sezhi, &
       &                          sigcx, scxlo, scxhi, sigcz, sczlo, sczhi) &
       bind(c,name='warpx_damp_pml_f_2d')
    integer, dimension(2), intent(in) :: tndlo, tndhi, flo, fhi
    integer, intent(in), value :: sexlo, sexhi, sezlo, sezhi, &
         &                        scxlo, scxhi, sczlo, sczhi
    real(amrex_real), intent(inout) :: f ( flo(1): fhi(1), flo(2): fhi(2),3)
    real(amrex_real), intent(in) :: sigex(sexlo:sexhi)
    real(amrex_real), intent(in) :: sigez(sezlo:sezhi)
    real(amrex_real), intent(in) :: sigcx(scxlo:scxhi)
    real(amrex_real), intent(in) :: sigcz(sczlo:sczhi)

    integer :: i,k

    do    k = tndlo(2), tndhi(2)
       do i = tndlo(1), tndhi(1)
          f(i,k,1) = f(i,k,1) * sigex(i)
          f(i,k,3) = f(i,k,3) * sigez(k)
       end do
    end do
  end subroutine warpx_damp_pml_f_2d


  subroutine warpx_damp_pml_f_3d (tndlo, tndhi, f, flo, fhi,&
       &                          sigex, sexlo, sexhi, sigey, seylo, seyhi, sigez, sezlo, sezhi, &
       &                          sigcx, scxlo, scxhi, sigcy, scylo, scyhi, sigcz, sczlo, sczhi) &
       bind(c,name='warpx_damp_pml_f_3d')
    integer, dimension(3), intent(in) :: tndlo, tndhi, flo, fhi
    integer, intent(in), value :: sexlo, sexhi, seylo, seyhi, sezlo, sezhi, &
         &                        scxlo, scxhi, scylo, scyhi, sczlo, sczhi
    real(amrex_real), intent(inout) :: f ( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),3)
    real(amrex_real), intent(in) :: sigex(sexlo:sexhi)
    real(amrex_real), intent(in) :: sigey(seylo:seyhi)
    real(amrex_real), intent(in) :: sigez(sezlo:sezhi)
    real(amrex_real), intent(in) :: sigcx(scxlo:scxhi)
    real(amrex_real), intent(in) :: sigcy(scylo:scyhi)
    real(amrex_real), intent(in) :: sigcz(sczlo:sczhi)

    integer :: i,j,k

    do       k = tndlo(3), tndhi(3)
       do    j = tndlo(2), tndhi(2)
          do i = tndlo(1), tndhi(1)
             f(i,j,k,1) = f(i,j,k,1) * sigex(i)
             f(i,j,k,2) = f(i,j,k,2) * sigey(j)
             f(i,j,k,3) = f(i,j,k,3) * sigez(k)
          end do
       end do
    end do
  end subroutine warpx_damp_pml_f_3d


end module warpx_pml_module
