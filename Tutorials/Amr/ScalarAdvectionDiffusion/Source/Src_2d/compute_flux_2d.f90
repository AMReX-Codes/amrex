module compute_flux_module

  implicit none

  private

  public :: &
       godunov_flux_2d, &
       mol2ndord_flux_2d, &
       mol4thord_flux_2d_limited, &
       mol4thord_flux_2d_nolimit, &
       getmaxmin_bounds

contains

  subroutine godunov_flux_2d(lo, hi, dt, dx, &
                             phi,ph_lo,ph_hi, &
                             umac,  u_lo,  u_hi, &
                             vmac,  v_lo,  v_hi, &
                             flxx, fx_lo, fx_hi, &
                             flxy, fy_lo, fy_hi, &
                             phix_1d, phiy_1d, phix, phiy, slope, glo, ghi, nu,  uselimit)

    use slope_module, only: slopex, slopey, slopex_nolim, slopey_nolim

    integer, intent(in) :: lo(2), hi(2), glo(2), ghi(2), uselimit
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

    if(uselimit .eq. 0) then
       call slopex_nolim(glo, ghi, &
            phi, ph_lo, ph_hi, &
            slope, glo, ghi)
    else

       call slopex(glo, ghi, &
            phi, ph_lo, ph_hi, &
            slope, glo, ghi)
    endif

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

    if(uselimit .eq. 0) then
       call slopey_nolim(glo, ghi, &
            phi, ph_lo, ph_hi, &
            slope, glo, ghi)
    else
       call slopey(glo, ghi, &
            phi, ph_lo, ph_hi, &
            slope, glo, ghi)
    endif
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
       phix_1d, phiy_1d, slope, glo, ghi, nu, uselimit)

    use slope_module, only: slopex, slopey, slopex_nolim, slopey_nolim

    integer, intent(in) :: lo(2), hi(2), glo(2), ghi(2), uselimit
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


    if(uselimit .eq. 0) then
       call slopex_nolim(glo, ghi, &
            phi, ph_lo, ph_hi, &
            slope, glo, ghi)
    else
       call slopex(glo, ghi, &
            phi, ph_lo, ph_hi, &
            slope, glo, ghi)
    endif

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
    
    if(uselimit .eq. 0) then
       call slopey_nolim(glo, ghi, &
            phi, ph_lo, ph_hi, &
            slope, glo, ghi)
    else
       call slopey(glo, ghi, &
            phi, ph_lo, ph_hi, &
            slope, glo, ghi)
    endif

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
  !velocity is pointwise on faces---no limiting here
  subroutine mol4thord_flux_2d_nolimit(lo, hi, dt, dx, &
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
         
    double precision :: diffflux, phicctemp, debtemp
    integer :: i, j
!    integer ::  numphi
!    double precision :: phitot

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
               ((27.0d0/24.0d0)*(phiptcc(i  ,j) - phiptcc(i-1,j)) &
               +( 1.0d0/24.0d0)*(phiptcc(i-2,j) - phiptcc(i+1,j)))

          fluxptx(i,j) = umac(i,j)*phiptx(i,j) + diffflux

!uncomment to just do diffusion
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
               ((27.0d0/24.0d0)*(phiptcc(i,j  ) - phiptcc(i,j-1)) &
               +( 1.0d0/24.0d0)*(phiptcc(i,j-2) - phiptcc(i,j+1)))

          fluxpty(i,j) = vmac(i,j)*phipty(i,j) + diffflux

!uncomment to just do diffusion
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

!debug   just set the flux = phiavex
!          flxx(i,j)  = phiavex(i,j)
!end debug
       end do
    end do

    do    j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
          debtemp  = fluxpty(i,j) + (1.0d0/24.d0)* &
               (fluxpty(i+1,j) + fluxpty(i-1,j) - 2.d0*fluxpty(i,j))
          flxy(i,j)  = debtemp

!debug   just set the flux = phiavex
!          flxy(i,j)  = phiavey(i,j)
!end debug
       end do
    end do

 
 end subroutine mol4thord_flux_2d_nolimit

  subroutine getmaxmin_bounds(lo, hi, glo, ghi, phi, maxval, minval, imax, imin)

    integer, intent(in) :: lo(2), hi(2), glo(2), ghi(2)
    double precision, intent(in   ) :: phi (glo(1):ghi(1),glo(2):ghi(2))
    double precision, intent(out   ) :: maxval, minval
    integer, intent(out) :: imax(2), imin(2)
    integer :: i, j

    maxval = -1.0d30
    minval =  1.0d30

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if(phi(i,j) .gt. maxval) then
             maxval = phi(i,j)
             imax(1) = i
             imax(2) = j
          endif
          if(phi(i,j) .lt. minval) then
             minval = phi(i,j)
             imin(1) = i
             imin(2) = j
          endif
       enddo
    enddo

    return
  end subroutine getmaxmin_bounds

  !phi coming in is assumed to be cell-averaged
  !velocity is pointwise on faces---this is the more complicated one that does limiting.
  subroutine mol4thord_flux_2d_limited(lo, hi, dt, dx, &
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

    use mempool_module, only : bl_allocate, bl_deallocate

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

    double precision, dimension(:,:), pointer, contiguous :: &
         dwfminux, dwfminuy, dwfplusx, dwfplusy, d2wfx, d2wfy, d2wcx, d2wcy, d3wx, d3wy, &
         wleftx, wrighx, wlefty, wrighy
         
         
    double precision :: diffflux, phicctemp, debtemp, mono1, mono2, signd2, d2wlim, &
         d2wc1, d2wc2, d2wc3, d2wf1, rho, d3wmin, d3wmax, dwthresh, maxw
    logical :: applylim
    integer :: i, j
    integer :: imax(2), imin(2), blo(2), bhi(2);
    double precision ::     maxval, minval
!    double precision :: phitot
!    integer ::  numphi


    call bl_allocate(dwfplusx,glo(1), ghi(1), glo(2), ghi(2))
    call bl_allocate(dwfminux,glo(1), ghi(1), glo(2), ghi(2))
    call bl_allocate(dwfplusy,glo(1), ghi(1), glo(2), ghi(2))
    call bl_allocate(dwfminuy,glo(1), ghi(1), glo(2), ghi(2))
    call bl_allocate(d2wfx   ,glo(1), ghi(1), glo(2), ghi(2))
    call bl_allocate(d2wcx   ,glo(1), ghi(1), glo(2), ghi(2))
    call bl_allocate(d3wx    ,glo(1), ghi(1), glo(2), ghi(2))
    call bl_allocate(d3wy    ,glo(1), ghi(1), glo(2), ghi(2))
    call bl_allocate(d2wfy   ,glo(1), ghi(1), glo(2), ghi(2))
    call bl_allocate(d2wcy   ,glo(1), ghi(1), glo(2), ghi(2))

    call bl_allocate(wleftx  ,glo(1), ghi(1), glo(2), ghi(2))
    call bl_allocate(wlefty  ,glo(1), ghi(1), glo(2), ghi(2))
    call bl_allocate(wrighx  ,glo(1), ghi(1), glo(2), ghi(2))
    call bl_allocate(wrighy  ,glo(1), ghi(1), glo(2), ghi(2))

    phiavey  = 1.0d30
    phiptcc  = 1.0d30
    dwfplusx = 1.0d30
    dwfminux = 1.0d30
    dwfplusy = 1.0d30
    dwfminuy = 1.0d30
    d2wfx    = 1.0d30
    d2wcx    = 1.0d30
    d3wx     = 1.0d30
    d2wfy    = 1.0d30
    d2wcy    = 1.0d30
    d3wy     = 1.0d30



    !STEP 0 
    ! 2.1 get cell-centered phi so we can compute a pointwise, fourth order gradient at faces
    ! needed  for diffusive fluxes (eqn 16 of m + c)


    do    j = lo(2)-3, hi(2)+3
       do i = lo(1)-3, hi(1)+3
          phicctemp  = phi(i,j) - (1.0d0/24.d0)* &
               (phi(i+1,j  )+phi(i-1,j  )&
               +phi(i  ,j+1)+phi(i  ,j-1)&
               -4.d0*phi(i,j))
          phiptcc(i,j)  = phicctemp
          
       end do
    end do

    blo(1) = lo(1)-3
    bhi(1) = hi(1)+3
    blo(2) = lo(2)-3
    bhi(2) = hi(2)+3
    call getmaxmin_bounds(blo, bhi, glo, ghi, phiptcc, maxval, minval, imax, imin)

    ! STEP ONE--- HYPERBOLIC FLUXES
    ! compute face average phi on x faces via eqn 17 of mccorquodale, colella
    
    do    j = lo(2)-2, hi(2)+2
       do i = lo(1)-1, hi(1)+2

          debtemp  = &
                (7.d0/12.d0)*(phi(i  ,j) + phi(i-1,j)) &
               -(1.d0/12.d0)*(phi(i-2,j) + phi(i+1,j))

          phiavex(i,j)  = debtemp

       end do
    end do

    blo(2) = lo(2)-2
    bhi(2) = hi(2)+2
    blo(1) = lo(1)-1
    bhi(1) = hi(1)+2
    call getmaxmin_bounds(blo, bhi, glo, ghi, phiavex, maxval, minval, imax, imin)
    
    !same for y faces
    do    j = lo(2)-1, hi(2)+2
       do i = lo(1)-2, hi(1)+2

          debtemp  = &
                (7.d0/12.d0)*(phi(i,j  ) + phi(i,j-1)) &
               -(1.d0/12.d0)*(phi(i,j-2) + phi(i,j+1))
          phiavey(i,j)  = debtemp

       end do
    end do

    blo(2) = lo(2)-1
    bhi(2) = hi(2)+2
    blo(1) = lo(1)-2
    bhi(1) = hi(1)+2
    call getmaxmin_bounds(blo, bhi, glo, ghi, phiavey, maxval, minval, imax, imin)

    !this stuff is for limiting 
    !The variable names are from mccorquodale & colella.   I could not come
    ! up with anything more sensible.   This is from section 2.4.1.
    do    j = lo(2)-2, hi(2)+2
       do i = lo(1)  , hi(1)
          dwfminux(i,j) =     phi(i  ,j) - phiavex(i,j)
          dwfplusx(i,j) = phiavex(i+1,j) -     phi(i,j)
       enddo
    enddo

    blo(2) = lo(2)-2
    bhi(2) = hi(2)+2
    blo(1) = lo(1)
    bhi(1) = hi(1)
    call getmaxmin_bounds(blo, bhi, glo, ghi, dwfminux, maxval, minval, imax, imin)
    call getmaxmin_bounds(blo, bhi, glo, ghi, dwfplusx, maxval, minval, imax, imin)

    do    j = lo(2)    , hi(2)
       do i = lo(1)-2  , hi(1)+2
          dwfminuy(i,j) =     phi(i,j  ) - phiavey(i,j)
          dwfplusy(i,j) = phiavey(i,j+1) -     phi(i,j)

       enddo
    enddo

    blo(2) = lo(2)
    bhi(2) = hi(2)
    blo(1) = lo(1)-2
    bhi(1) = hi(1)+2
    call getmaxmin_bounds(blo, bhi, glo, ghi, dwfminuy, maxval, minval, imax, imin)
    call getmaxmin_bounds(blo, bhi, glo, ghi, dwfplusy, maxval, minval, imax, imin)

    do    j = lo(2)-2, hi(2)+2
       do i = lo(1)-1, hi(1)+1
          d2wfx(i,j) = 6.0d0*(phiavex(i  ,j) + phiavex(i+1,j) - 2.0d0*phi(i,j))
          d2wcx(i,j) =       (    phi(i-1,j) +     phi(i+1,j) - 2.0d0*phi(i,j))

       enddo
    enddo


    blo(2) = lo(2)-2
    bhi(2) = hi(2)+2
    blo(1) = lo(1)-1
    bhi(1) = hi(1)+1
    call getmaxmin_bounds(blo, bhi, glo, ghi, d2wfx, maxval, minval, imax, imin)
    call getmaxmin_bounds(blo, bhi, glo, ghi, d2wcx, maxval, minval, imax, imin)

    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)-2, hi(1)+2
          d2wfy(i,j) = 6.0d0*(phiavey(i,j  ) + phiavey(i,j+1) - 2.0d0*phi(i,j))
          d2wcy(i,j) =       (    phi(i,j-1) +     phi(i,j+1) - 2.0d0*phi(i,j))

       enddo
    enddo

    blo(2) = lo(2)-1
    bhi(2) = hi(2)+1
    blo(1) = lo(1)-2
    bhi(1) = hi(1)+2
    call getmaxmin_bounds(blo, bhi, glo, ghi, d2wfy, maxval, minval, imax, imin)
    call getmaxmin_bounds(blo, bhi, glo, ghi, d2wcy, maxval, minval, imax, imin)


    do    j = lo(2)-2, hi(2)+2
       do i = lo(1),   hi(1)+1
          d3wx(i,j) = d2wcx(i,j) - d2wcx(i-1,j)
       enddo
    enddo


    blo(2) = lo(2)-2
    bhi(2) = hi(2)+2
    blo(1) = lo(1)
    bhi(1) = hi(1)+1
    call getmaxmin_bounds(blo, bhi, glo, ghi, d3wx, maxval, minval, imax, imin)


    do    j = lo(2)  , hi(2)+1
       do i = lo(1)-2, hi(1)+2
          d3wy(i,j) = d2wcy(i,j) - d2wcy(i,j-1)
       enddo
    enddo


    blo(2) = lo(2)
    bhi(2) = hi(2)+1
    blo(1) = lo(1)-2
    bhi(1) = hi(1)+2
    call getmaxmin_bounds(blo, bhi, glo, ghi, d3wy, maxval, minval, imax, imin)


    !initialize left and right states
    do    j = lo(2)-2, hi(2)+2
       do i = lo(1)  , hi(1)+1
          wleftx(i,j) = phiavex(i,j)
          wrighx(i,j) = phiavex(i,j)

       enddo
    enddo
    do    j = lo(2)  , hi(2)+1
       do i = lo(1)-2, hi(1)+2
          wlefty(i,j) = phiavey(i,j)
          wrighy(i,j) = phiavey(i,j)
       enddo
    enddo

    !now for the really complicted stuff in 2.4.1
    !If there are extrema detected a cell, modify the states
    !that are extrapolate from that cell.  If monotonoic but there
    ! is a huge jump in the state, mollify that
    do    j = lo(2)-2, hi(2)+2
       do i = lo(1)  , hi(1)+1
          mono1 = dwfminux(i,j) * dwfplusx(i,j)
          mono2 = (phi(i,j)-phi(i-2,j))*(phi(i+2,j)-phi(i,j))
          if((mono1 .lt. 0.0d0).or.(mono2 .lt. 0.0d0)) then

             !see if all second derivs in sight are the same sign.   If they are,
             ! calculate d2wlim.     Otherwise, d2wlim = 0 
             !this stuff comes from equation 26
             d2wlim = 0.0d0
             signd2 = 0.0d0
             d2wc1 = d2wcx(i-1,j)
             d2wc2 = d2wcx(i  ,j)
             d2wc3 = d2wcx(i+1,j)
             d2wf1 = d2wfx(i  ,j)
             d2wlim = min(abs(d2wf1), 1.25d0*abs(d2wc1))
             d2wlim = min(d2wlim    , 1.25d0*abs(d2wc2))
             d2wlim = min(d2wlim    , 1.25d0*abs(d2wc3))
             if((d2wc1.gt.0.0d0) .and. &
                (d2wc2.gt.0.0d0) .and. &
                (d2wc3.gt.0.0d0) .and. &
                (d2wf1.gt.0.0d0)) then
                signd2 = 1.0d0
                d2wlim = d2wlim*signd2
             else if((d2wc1.lt.0.0d0) .and. &
                     (d2wc2.lt.0.0d0) .and. &
                     (d2wc3.lt.0.0d0) .and. &
                     (d2wf1.lt.0.0d0)) then
                signd2 = -1.0d0
                d2wlim = d2wlim*signd2
             else
                d2wlim = 0.0d0
             endif
             !equation 27 below
             maxw = max(abs(phi(i,j)),abs(phi(i+1,j)))
             maxw = max(maxw         ,abs(phi(i-1,j)))
             maxw = max(maxw         ,abs(phi(i+2,j)))
             maxw = max(maxw         ,abs(phi(i-2,j)))
             dwthresh = 1.0d-12*(maxw)
             if(abs(d2wf1).lt.dwthresh) then
                rho = 0.0d0
             else
                rho = d2wlim/d2wf1
             endif
             if(rho .lt. (1.0d0 - 1.0d-12)) then
                applylim = .false.
             else
                !this stuff is equation 28
                d3wmin = min(d3wx(i,j),d3wx(i-1,j))
                d3wmin = min(d3wmin   ,d3wx(i+1,j))
                d3wmin = min(d3wmin   ,d3wx(i+2,j))
                d3wmax = max(d3wx(i,j),d3wx(i-1,j))
                d3wmax = max(d3wmax   ,d3wx(i+1,j))
                d3wmax = max(d3wmax   ,d3wx(i+2,j))
                dwthresh = 0.1d0*max(abs(d3wmax),abs(d3wmin))
                if((d3wmax - d3wmin).lt. dwthresh) then
                   applylim = .false.
                else
                   applylim = .true.
                endif
             endif
             !the following covered equations 29-32
             if(applylim) then
                if((dwfminux(i,j)*dwfplusx(i,j)).lt.0.0) then

                   wrighx(i,j) = phi(i,j) -  rho*dwfminux(i,j) !eqn 29 (corrected)
                   wleftx(i,j) = phi(i,j) +  rho*dwfplusx(i,j) !eqn 30

                else if(abs(dwfminux(i,j)) .ge. 2.0d0*abs(dwfplusx(i,j))) then
                   !equation 31
                   wrighx(i,j) = phi(i,j) -  2.0d0*(1.0d0-rho)*dwfplusx(i,j) - rho*dwfminux(i,j)
                else if(abs(dwfplusx(i,j)) .ge. 2.0d0*abs(dwfminux(i,j))) then
                   !equation 32
                   wleftx(i,j) = phi(i,j) +  2.0d0*(1.0d0-rho)*dwfminux(i,j) + rho*dwfplusx(i,j)
                endif
             endif
          else !monotonic but still need to check for big jumps in the state
             if(abs(dwfminux(i,j)) .ge. (2.0d0*abs(dwfplusx(i,j)))) then
                wrighx(i,j) = phi(i,j) -  2.0d0*dwfplusx(i,j) !eqn 33
             endif
             if(abs(dwfplusx(i,j)) .ge. (2.0d0*abs(dwfminux(i,j)))) then
                wleftx(i,j) = phi(i,j) +  2.0d0*dwfminux(i,j) !eqn 34
             endif
          endif
       enddo
    enddo
    !now we do it all again for the y direction. 
    do    j = lo(2)  , hi(2)+1
       do i = lo(1)-2, hi(1)+2
          mono1 = dwfminuy(i,j) * dwfplusy(i,j)
          mono2 = (phi(i,j)-phi(i,j-2))*(phi(i,j+2)-phi(i,j))
          if((mono1 .lt. 0.0d0).or.(mono2 .lt. 0.0d0)) then

             !see if all second derivs in sight are the same sign.   If they are,
             ! calculate d2wlim.     Otherwise, d2wlim = 0 
             !this stuff comes from equation 26
             d2wlim = 0.0d0
             signd2 = 0.0d0
             d2wc1 = d2wcy(i,j-1)
             d2wc2 = d2wcy(i,j  )
             d2wc3 = d2wcy(i,j+1)
             d2wf1 = d2wfy(i,j  )
             d2wlim = min(abs(d2wf1), 1.25d0*abs(d2wc1))
             d2wlim = min(d2wlim    , 1.25d0*abs(d2wc2))
             d2wlim = min(d2wlim    , 1.25d0*abs(d2wc3))
             if((d2wc1.gt.0.0d0) .and. &
                (d2wc2.gt.0.0d0) .and. &
                (d2wc3.gt.0.0d0) .and. &
                (d2wf1.gt.0.0d0)) then
                signd2 = 1.0d0
                d2wlim = d2wlim*signd2
             else if((d2wc1.lt.0.0d0) .and. &
                     (d2wc2.lt.0.0d0) .and. &
                     (d2wc3.lt.0.0d0) .and. &
                     (d2wf1.lt.0.0d0)) then
                signd2 = -1.0d0
                d2wlim = d2wlim*signd2
             else
                d2wlim = 0.0d0
             endif
             !equation 27 below
             maxw = max(abs(phi(i,j)),abs(phi(i,j+1)))
             maxw = max(maxw         ,abs(phi(i,j-1)))
             maxw = max(maxw         ,abs(phi(i,j+2)))
             maxw = max(maxw         ,abs(phi(i,j-2)))
             dwthresh = 1.0d-12*(maxw)
             if(abs(d2wf1).lt.dwthresh) then
                rho = 0.0d0
             else
                rho = d2wlim/d2wf1
             endif
             if(rho .lt. (1.0d0 - 1.0d-12)) then
                applylim = .false.
             else
                !this stuff is equation 28
                d3wmin = min(d3wy(i,j),d3wy(i,j-1))
                d3wmin = min(d3wmin   ,d3wy(i,j+1))
                d3wmin = min(d3wmin   ,d3wy(i,j+2))
                d3wmax = max(d3wy(i,j),d3wy(i,j-1))
                d3wmax = max(d3wmax   ,d3wy(i,j+1))
                d3wmax = max(d3wmax   ,d3wy(i,j+2))
                dwthresh = 0.1d0*max(abs(d3wmax),abs(d3wmin))
                if((d3wmax - d3wmin).lt. dwthresh) then
                   applylim = .false.
                else
                   applylim = .true.
                endif
             endif
             !the following covered equations 29-32
             if(applylim) then
                if((dwfminuy(i,j)*dwfplusy(i,j)).lt.0.0) then

                   wrighy(i,j) = phi(i,j) -  rho*dwfminuy(i,j) !eqn 29 (corrected)
                   wlefty(i,j) = phi(i,j) +  rho*dwfplusy(i,j) !eqn 30

                else if(abs(dwfminuy(i,j)) .ge. 2.0d0*abs(dwfplusy(i,j))) then
                   !equation 31
                   wrighy(i,j) = phi(i,j) -  2.0d0*(1.0d0-rho)*dwfplusy(i,j) - rho*dwfminuy(i,j)
                else if(abs(dwfplusy(i,j)) .ge. 2.0d0*abs(dwfminuy(i,j))) then
                   !equation 32
                   wlefty(i,j) = phi(i,j) +  2.0d0*(1.0d0-rho)*dwfminuy(i,j) + rho*dwfplusy(i,j)
                endif
             endif

          else !monotonic but still need to check for big jumps in the state
             if(abs(dwfminuy(i,j)) .ge. (2.0d0*abs(dwfplusy(i,j)))) then
                wrighy(i,j) = phi(i,j) -  2.0d0*dwfplusy(i,j) !eqn 33
             endif
             if(abs(dwfplusx(i,j)) .ge. (2.0d0*abs(dwfminux(i,j)))) then
                wlefty(i,j) = phi(i,j) +  2.0d0*dwfminuy(i,j) !eqn 34
             endif
          endif
    
       enddo
    enddo
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
               ((27.0d0/24.0d0)*(phiptcc(i  ,j) - phiptcc(i-1,j)) &
               +( 1.0d0/24.0d0)*(phiptcc(i-2,j) - phiptcc(i+1,j)))

          fluxptx(i,j) = umac(i,j)*phiptx(i,j) + diffflux

!uncomment to just do diffusion
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
               ((27.0d0/24.0d0)*(phiptcc(i,j  ) - phiptcc(i,j-1)) &
               +( 1.0d0/24.0d0)*(phiptcc(i,j-2) - phiptcc(i,j+1)))

          fluxpty(i,j) = vmac(i,j)*phipty(i,j) + diffflux

!uncomment to just do diffusion
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

!debug   just set the flux = phiavex
!          flxx(i,j)  = phiavex(i,j)
!end debug
       end do
    end do

    do    j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
          debtemp  = fluxpty(i,j) + (1.0d0/24.d0)* &
               (fluxpty(i+1,j) + fluxpty(i-1,j) - 2.d0*fluxpty(i,j))
          flxy(i,j)  = debtemp

!debug   just set the flux = phiavex
!          flxy(i,j)  = phiavey(i,j)
!end debug
       end do
    end do

    call bl_deallocate(dwfplusx)
    call bl_deallocate(dwfminux)
    call bl_deallocate(dwfplusy)
    call bl_deallocate(dwfminuy)
    call bl_deallocate(d2wfx   )
    call bl_deallocate(d2wcx   )
    call bl_deallocate(d3wx    )
    call bl_deallocate(d3wy    )
    call bl_deallocate(d2wfy   )
    call bl_deallocate(d2wcy   )
    call bl_deallocate(d3wy    )
    call bl_deallocate(wleftx  )
    call bl_deallocate(wlefty  )
    call bl_deallocate(wrighx  )
    call bl_deallocate(wrighy  )

 
 end subroutine mol4thord_flux_2d_limited

end module compute_flux_module
