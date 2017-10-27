module compute_flux_module

  implicit none

  private

  public :: godunov_flux_2d, mol2ndord_flux_2d, mol4thord_flux_2d

contains

  subroutine godunov_flux_2d(lo, hi, dt, dx, &
                             phi,ph_lo,ph_hi, &
                             umac,  u_lo,  u_hi, &
                             vmac,  v_lo,  v_hi, &
                             flxx, fx_lo, fx_hi, &
                             flxy, fy_lo, fy_hi, &
                             phix_1d, phiy_1d, phix, phiy, slope, glo, ghi, nu)

    use slope_module, only: slopex, slopey

    integer, intent(in) :: lo(2), hi(2), glo(2), ghi(2)
    double precision, intent(in) :: dt, dx(2), nu
    integer, intent(in) :: ph_lo(2), ph_hi(2)
    integer, intent(in) ::  u_lo(2),  u_hi(2)
    integer, intent(in) ::  v_lo(2),  v_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    double precision, intent(in   ) :: phi (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2))
    double precision, intent(in   ) :: umac( u_lo(1): u_hi(1), u_lo(2): u_hi(2))
    double precision, intent(in   ) :: vmac( v_lo(1): v_hi(1), v_lo(2): v_hi(2))
    double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
    double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))
    double precision, dimension(glo(1):ghi(1),glo(2):ghi(2)) :: &
         phix_1d, phiy_1d, phix, phiy, slope
         
    integer :: i, j
    double precision :: hdtdx(2)

    hdtdx = 0.5*(dt/dx)

    call slopex(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)

    ! compute phi on x faces using umac to upwind; ignore transverse terms
    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)  , hi(1)+1

          if (umac(i,j) .lt. 0.d0) then
             phix_1d(i,j) = phi(i  ,j) - (0.5d0 + hdtdx(1)*umac(i,j))*slope(i  ,j)
          else
             phix_1d(i,j) = phi(i-1,j) + (0.5d0 - hdtdx(1)*umac(i,j))*slope(i-1,j)
          end if

       end do
    end do

    call slopey(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)

    ! compute phi on y faces using umac to upwind; ignore transverse terms
    do    j = lo(2)  , hi(2)+1
       do i = lo(1)-1, hi(1)+1

          if (vmac(i,j) .lt. 0.d0) then
             phiy_1d(i,j) = phi(i,j  ) - (0.5d0 + hdtdx(2)*vmac(i,j))*slope(i,j  )
          else
             phiy_1d(i,j) = phi(i,j-1) + (0.5d0 - hdtdx(2)*vmac(i,j))*slope(i,j-1)
          end if

       end do
    end do

    ! update phi on x faces by adding in y-transverse terms
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)+1

          if (umac(i,j) .lt. 0.d0) then
             phix(i,j) = phix_1d(i,j) &
                  - hdtdx(2)*( 0.5d0*(vmac(i  ,j+1)+vmac(i  ,j)) * (phiy_1d(i  ,j+1)-phiy_1d(i  ,j)) )
          else
             phix(i,j) = phix_1d(i,j) &
                  - hdtdx(2)*( 0.5d0*(vmac(i-1,j+1)+vmac(i-1,j)) * (phiy_1d(i-1,j+1)-phiy_1d(i-1,j)) )
          end if

          ! compute final x-fluxes
          ! including diffusive fluxes
          flxx(i,j) = phix(i,j)*umac(i,j)  - nu*(phi(i,j) - phi(i-1,j))/dx(1)

       end do
    end do

    ! update phi on y faces by adding in x-transverse terms
    do    j = lo(2), hi(2)+1
       do i = lo(1), hi(1)

          if (vmac(i,j) .lt. 0.d0) then
             phiy(i,j) = phiy_1d(i,j) &
                  - hdtdx(1)*( 0.5d0*(umac(i+1,j  )+umac(i,j  )) * (phix_1d(i+1,j  )-phix_1d(i,j  )) )
          else
             phiy(i,j) = phiy_1d(i,j) &
                  - hdtdx(1)*( 0.5d0*(umac(i+1,j-1)+umac(i,j-1)) * (phix_1d(i+1,j-1)-phix_1d(i,j-1)) )
          end if

          ! compute final y-fluxes
          ! including diffusive fluxes
          flxy(i,j) = phiy(i,j)*vmac(i,j)  - nu*(phi(i,j) - phi(i,j-1))/dx(2)

       end do
    end do

  end subroutine godunov_flux_2d


  subroutine mol2ndord_flux_2d(lo, hi, dt, dx, &
       phi,ph_lo,ph_hi, &
       umac,  u_lo,  u_hi, &
       vmac,  v_lo,  v_hi, &
       flxx, fx_lo, fx_hi, &
       flxy, fy_lo, fy_hi, &
       phix_1d, phiy_1d, slope, glo, ghi, nu)

    use slope_module, only: slopex, slopey

    integer, intent(in) :: lo(2), hi(2), glo(2), ghi(2)
    double precision, intent(in) :: dt, dx(2), nu
    integer, intent(in) :: ph_lo(2), ph_hi(2)
    integer, intent(in) ::  u_lo(2),  u_hi(2)
    integer, intent(in) ::  v_lo(2),  v_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    double precision, intent(in   ) :: phi (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2))
    double precision, intent(in   ) :: umac( u_lo(1): u_hi(1), u_lo(2): u_hi(2))
    double precision, intent(in   ) :: vmac( v_lo(1): v_hi(1), v_lo(2): v_hi(2))
    double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
    double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))
    double precision, dimension(glo(1):ghi(1),glo(2):ghi(2)) :: &
         phix_1d, phiy_1d, slope
    double precision phi1, phi2, ftemp, veltemp, diffterm
         
    integer :: i, j


    call slopex(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)

    ! compute phi on x faces using umac to upwind; ignore transverse terms
    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)  , hi(1)+1

          phi1 =  phi(i  ,j) - 0.5d0*slope(i  ,j)
          phi2 =  phi(i-1,j) + 0.5d0*slope(i-1,j)

          if (umac(i,j) .lt. 0.d0) then
             phix_1d(i,j) = phi(i  ,j) - 0.5d0*slope(i  ,j)
          else
             phix_1d(i,j) = phi(i-1,j) + 0.5d0*slope(i-1,j)
          end if

       end do
    end do

    call slopey(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)

    ! compute phi on y faces using umac to upwind; ignore transverse terms
    do    j = lo(2)  , hi(2)+1
       do i = lo(1)-1, hi(1)+1

          if (vmac(i,j) .lt. 0.d0) then
             phiy_1d(i,j) = phi(i,j  ) - (0.5d0)*slope(i,j  )
          else
             phiy_1d(i,j) = phi(i,j-1) + (0.5d0)*slope(i,j-1)
          end if

       end do
    end do

    ! compute final x-fluxes
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)+1
          ! including diffusive fluxes
          veltemp  = umac(i,j) 
          diffterm = - nu*(phi(i,j) - phi(i-1,j))/dx(1)
          phi1  = phix_1d(i,j)
          ftemp  = phix_1d(i,j)*umac(i,j)  - nu*(phi(i,j) - phi(i-1,j))/dx(1)
          flxx(i,j) = phix_1d(i,j)*umac(i,j)  - nu*(phi(i,j) - phi(i-1,j))/dx(1)
       end do
    end do

    ! compute final y-fluxes
    do    j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
          ! including diffusive fluxes
          ftemp  = phiy_1d(i,j)*vmac(i,j)  - nu*(phi(i,j) - phi(i,j-1))/dx(2)

          flxy(i,j) = phiy_1d(i,j)*vmac(i,j)  - nu*(phi(i,j) - phi(i,j-1))/dx(2)
       end do
    end do

  end subroutine mol2ndord_flux_2d



  !phi coming in is assumed to be cell-averaged
  !velocity is pointwise on faces
  subroutine mol4thord_flux_2d(lo, hi, dt, dx, &
                               phi ,ph_lo,ph_hi, &
                               umac,  u_lo,  u_hi, &
                               vmac,  v_lo,  v_hi, &
                               flxx, fx_lo, fx_hi, &
                               flxy, fy_lo, fy_hi, &
                               fluxptx, phiptx, phiavex, &
                               fluxpty, phipty, phiavey, &
                               phiptcc, glo, ghi, nu, &
                               deblocell, debhicell, &
                               hisidedebfacelo, hisidedebfacehi, &
                               losidedebfacelo, losidedebfacehi, printstuff)

    integer, intent(in) :: lo(2), hi(2), glo(2), ghi(2), printstuff
    integer, intent(in) :: deblocell(2), debhicell(2)
    integer, intent(in) :: hisidedebfacelo(2), hisidedebfacehi(2)
    integer, intent(in) :: losidedebfacelo(2), losidedebfacehi(2)
    double precision, intent(in) :: dt, dx(2), nu
    integer, intent(in) :: ph_lo(2), ph_hi(2)
    integer, intent(in) ::  u_lo(2),  u_hi(2)
    integer, intent(in) ::  v_lo(2),  v_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    double precision, intent(in   ) :: phi (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2))
    double precision, intent(in   ) :: umac( u_lo(1): u_hi(1), u_lo(2): u_hi(2))
    double precision, intent(in   ) :: vmac( v_lo(1): v_hi(1), v_lo(2): v_hi(2))
    double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
    double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))
    double precision, dimension(glo(1):ghi(1),glo(2):ghi(2)) :: &
         fluxptx, fluxpty,  phiptx, phipty,  phiavex, phiavey, phiptcc
         
    double precision :: diffflux, phicctemp, debtemp, phitot
    integer :: i, j, numphi

    fluxptx = 1.0d30
    fluxpty = 1.0d30
    phiptx  = 1.0d30
    phipty  = 1.0d30
    phiavex = 1.0d30
    phiavey = 1.0d30
    phiptcc = 1.0d30
    !STEP 0 
    ! 2.1 get cell-centered phi so we can compute a pointwise, fourth order gradient at faces
    ! needed  for diffusive fluxes

    do    j = lo(2)-3, hi(2)+3
       do i = lo(1)-3, hi(1)+3
          phicctemp  = phi(i,j) - (1.0d0/24.d0)* &
               (phi(i+1,j  )+phi(i-1,j  )&
               +phi(i  ,j+1)+phi(i  ,j-1)&
               -4.d0*phi(i,j))
          phiptcc(i,j)  = phicctemp
       end do
    end do

    ! STEP ONE--- HYPERBOLIC FLUXES
    ! compute face average phi on x faces via eqn 17 of mccorquodale, colella
    
    do    j = lo(2)-2, hi(2)+2
       do i = lo(1)  , hi(1)+1

          debtemp  = &
                (7.d0/12.d0)*(phi(i  ,j) + phi(i-1,j)) &
               -(1.d0/12.d0)*(phi(i-2,j) + phi(i+1,j))

          phiavex(i,j)  = debtemp
       end do
    end do

    
    !same for y faces
    do    j = lo(2)  , hi(2)+1
       do i = lo(1)-2, hi(1)+2

          debtemp  = &
                (7.d0/12.d0)*(phi(i,j  ) + phi(i,j-1)) &
               -(1.d0/12.d0)*(phi(i,j-2) + phi(i,j+1))
          phiavey(i,j)  = debtemp

       end do
    end do

    !now get point valued phi at faces so we can multiply by point valued velocity
    ! phipt = phiave - (h^2/24)*(lapl^2d(phi_ave))
    ! while I am at it, multiply in pointwise velocity so we get pointwise *HYPERBOLIC* flux
    !also  get pointwise diffusive fluxes using 4th order finite differences of pointwise phi

    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)  , hi(1)+1
          phiptx(i,j)  =    phiavex(i,j  ) - (1.0d0/24.d0)* &
               (phiavex(i,j+1) + phiavex(i,j-1) - 2.d0*phiavex(i,j))

          debtemp  = phiptx(i,j)

          diffflux  = (-nu/dx(1))* &
               ( phiptcc(i  ,j) - phiptcc(i-1,j) &
               -(phiptcc(i+1,j) + phiptcc(i-1,j) - 2.0d0*phiptcc(i  ,j)) &
               +(phiptcc(i  ,j) + phiptcc(i-2,j) - 2.0d0*phiptcc(i-1,j)))

          fluxptx(i,j) = umac(i,j)*phiptx(i,j) + diffflux
!          fluxptx(i,j) =  diffflux

       end do
    end do

    
    !same for y faces
    do    j = lo(2)  , hi(2)+1
       do i = lo(1)-1, hi(1)+1

          phipty(i,j)  = phiavey(i,j) - (1.0d0/24.d0)* &
               (phiavey(i+1,j)+phiavey(i-1,j)-2.d0*phiavey(i,j))

          debtemp  = phipty(i,j)

          diffflux  = (-nu/dx(2))* &
               ( phiptcc(i,j  ) - phiptcc(i,j-1) &
               -(phiptcc(i,j+1) + phiptcc(i,j-1) - 2.0d0*phiptcc(i,j  )) &
               +(phiptcc(i,j  ) + phiptcc(i,j-2) - 2.0d0*phiptcc(i,j-1)))

          fluxpty(i,j) = vmac(i,j)*phipty(i,j) + diffflux
!          fluxpty(i,j) =  diffflux

       end do
    end do

    ! now transform pointwise fluxes into face-averaged fluxes
    !fluxave = fluxpt + (1/24)(Lapl2d(fluxpt))
    
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)+1
          debtemp  = fluxptx(i,j) + (1.0d0/24.d0)* &
               (fluxptx(i,j+1) + fluxptx(i,j-1) - 2.d0*fluxptx(i,j))

          flxx(i,j)  = debtemp
       end do
    end do

    do    j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
          debtemp  = fluxpty(i,j) + (1.0d0/24.d0)* &
               (fluxpty(i+1,j) + fluxpty(i-1,j) - 2.d0*fluxpty(i,j))
          flxy(i,j)  = debtemp

       end do
    end do

 
 end subroutine mol4thord_flux_2d


end module compute_flux_module
