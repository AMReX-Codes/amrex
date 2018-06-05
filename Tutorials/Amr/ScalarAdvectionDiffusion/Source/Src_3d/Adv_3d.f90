
subroutine advectDiffGodunov(time, lo, hi, &
     &            uin , ui_lo, ui_hi, &
     &            dphdtout, uo_lo, uo_hi, &
     &            vx  , vx_lo, vx_hi, &
     &            vy  , vy_lo, vy_hi, &
     &            vz  , vz_lo, vz_hi, &
     &            flxx, fx_lo, fx_hi, &
     &            flxy, fy_lo, fy_hi, &
     &            flxz, fz_lo, fz_hi, &
     &            dx,dt, nu, uselimit) bind(C, name="advectDiffGodunov")
  
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use compute_flux_module, only : godunov_flux_3d

  implicit none

  integer, intent(in) :: lo(3), hi(3),  uselimit
  double precision, intent(in) :: dx(3), dt, time, nu
  integer, intent(in) :: ui_lo(3), ui_hi(3)
  integer, intent(in) :: uo_lo(3), uo_hi(3)
  integer, intent(in) :: vx_lo(3), vx_hi(3)
  integer, intent(in) :: vy_lo(3), vy_hi(3)
  integer, intent(in) :: vz_lo(3), vz_hi(3)
  integer, intent(in) :: fx_lo(3), fx_hi(3)
  integer, intent(in) :: fy_lo(3), fy_hi(3)
  integer, intent(in) :: fz_lo(3), fz_hi(3)
  double precision, intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3))
  double precision, intent(inout) :: dphdtout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3))
  double precision, intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2),vx_lo(3):vx_hi(3))
  double precision, intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2),vy_lo(3):vy_hi(3))
  double precision, intent(in   ) :: vz  (vz_lo(1):vz_hi(1),vz_lo(2):vz_hi(2),vz_lo(3):vz_hi(3))
  double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
  double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
  double precision, intent(  out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))

  integer :: i, j, k
  integer :: glo(3), ghi(3)
  double precision ::  umax, vmax, wmax

  ! Some compiler may not support 'contiguous'.  Remove it in that case.
  double precision, dimension(:,:,:), pointer, contiguous :: &
       phix, phix_y, phix_z, phiy, phiy_x, phiy_z, phiz, phiz_x, phiz_y, slope


  glo = lo - 1
  ghi = hi + 1

  ! edge states
  call bl_allocate(phix  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phix_y,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phix_z,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiy  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiy_x,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiy_z,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiz  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiz_x,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiz_y,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  ! slope
  call bl_allocate(slope,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))  
  
  ! We like to allocate these **pointers** here and then pass them to a function
  ! to remove their pointerness for performance, because normally pointers could
  ! be aliasing.  We need to use pointers instead of allocatable arrays because
  ! we like to use BoxLib's bl_allocate to allocate memeory instead of the intrinsic
  ! allocate.  Bl_allocate is much faster than allocate inside OMP.  
  ! Note that one MUST CALL BL_DEALLOCATE.

  ! check if CFL condition is violated.
  umax = maxval(abs(vx))
  vmax = maxval(abs(vy))
  wmax = maxval(abs(vz))
  if ( umax*dt .ge. dx(1) .or. &
       vmax*dt .ge. dx(2) .or. &
       wmax*dt .ge. dx(3) ) then
     print *, "umax = ", umax, ", vmax = ", vmax, ", wmax = ", wmax, ", dt = ", dt, ", dx = ", dx
     call bl_error("CFL violation. Use smaller adv.cfl.")
  end if

  ! call a function to compute flux
  call godunov_flux_3d(lo, hi, dt, dx, &
                       uin, ui_lo, ui_hi, &
                       vx, vx_lo, vx_hi, &
                       vy, vy_lo, vy_hi, &
                       vz, vz_lo, vz_hi, &
                       flxx, fx_lo, fx_hi, &
                       flxy, fy_lo, fy_hi, &
                       flxz, fz_lo, fz_hi, &
                       phix, phix_y, phix_z, &
                       phiy, phiy_x, phiy_z, &
                       phiz, phiz_x, phiz_y, &
                       slope, glo, ghi, nu, uselimit)

  ! Do a conservative update
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           !notice that some mad scientist has reversed the flux difference
           ! rather than just put a negative sign out front.
           dphdtout(i,j,k) =  &
                ( (flxx(i,j,k) - flxx(i+1,j  ,k  )) /dx(1) &
                + (flxy(i,j,k) - flxy(i  ,j+1,k  )) /dx(2) &
                + (flxz(i,j,k) - flxz(i  ,j,  k+1)) /dx(3) )
        enddo
     enddo
  enddo

  ! Scale by face area in order to correctly reflux because flux register requires this
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
           flxx(i,j,k) = flxx(i,j,k) * (dt * dx(2)*dx(3))
        enddo
     enddo
  enddo
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)+1 
        do i = lo(1), hi(1)
           flxy(i,j,k) = flxy(i,j,k) * (dt * dx(1)*dx(3))
        enddo
     enddo
  enddo
  do       k = lo(3), hi(3)+1
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           flxz(i,j,k) = flxz(i,j,k) * (dt * dx(1)*dx(2))
        enddo
     enddo
  enddo

  call bl_deallocate(phix  )
  call bl_deallocate(phix_y)
  call bl_deallocate(phix_z)
  call bl_deallocate(phiy  )
  call bl_deallocate(phiy_x)
  call bl_deallocate(phiy_z)
  call bl_deallocate(phiz  )
  call bl_deallocate(phiz_x)
  call bl_deallocate(phiz_y)
  call bl_deallocate(slope)

end subroutine advectDiffGodunov



subroutine advectDiffMOL2ndOrd(time, lo, hi, &
     &            uin , ui_lo, ui_hi, &
     &            dphdtout, uo_lo, uo_hi, &
     &            vx  , vx_lo, vx_hi, &
     &            vy  , vy_lo, vy_hi, &
     &            vz  , vz_lo, vz_hi, &
     &            flxx, fx_lo, fx_hi, &
     &            flxy, fy_lo, fy_hi, &
     &            flxz, fz_lo, fz_hi, &
     &            dx,dt,nu,uselimit) bind(C, name="advectDiffMOL2ndOrd")
  
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use compute_flux_module, only : mol2ndord_flux_3d
  implicit none

  integer, intent(in) :: lo(3), hi(3), uselimit
  double precision, intent(in) :: dx(3), dt, time, nu
  integer, intent(in) :: ui_lo(3), ui_hi(3)
  integer, intent(in) :: uo_lo(3), uo_hi(3)
  integer, intent(in) :: vx_lo(3), vx_hi(3)
  integer, intent(in) :: vy_lo(3), vy_hi(3)
  integer, intent(in) :: vz_lo(3), vz_hi(3)
  integer, intent(in) :: fx_lo(3), fx_hi(3)
  integer, intent(in) :: fy_lo(3), fy_hi(3)
  integer, intent(in) :: fz_lo(3), fz_hi(3)
  double precision, intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3))
  double precision, intent(inout) :: dphdtout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3))
  double precision, intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2),vx_lo(3):vx_hi(3))
  double precision, intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2),vy_lo(3):vy_hi(3))
  double precision, intent(in   ) :: vz  (vz_lo(1):vz_hi(1),vz_lo(2):vz_hi(2),vz_lo(3):vz_hi(3))
  double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
  double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
  double precision, intent(  out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))

  integer :: i, j, k
  integer :: glo(3), ghi(3)
  double precision ::  umax, vmax, wmax

  ! Some compiler may not support 'contiguous'.  Remove it in that case.
  double precision, dimension(:,:,:), pointer, contiguous :: &
       phix, phiy, phiz, slope


  glo = lo - 1
  ghi = hi + 1

  ! edge states
  call bl_allocate(phix  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiy  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiz  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  ! slope
  call bl_allocate(slope,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))  
  
  ! We like to allocate these **pointers** here and then pass them to a function
  ! to remove their pointerness for performance, because normally pointers could
  ! be aliasing.  We need to use pointers instead of allocatable arrays because
  ! we like to use BoxLib's bl_allocate to allocate memeory instead of the intrinsic
  ! allocate.  Bl_allocate is much faster than allocate inside OMP.  
  ! Note that one MUST CALL BL_DEALLOCATE.

  ! check if CFL condition is violated.
  umax = maxval(abs(vx))
  vmax = maxval(abs(vy))
  wmax = maxval(abs(vz))
  if ( umax*dt .ge. dx(1) .or. &
       vmax*dt .ge. dx(2) .or. &
       wmax*dt .ge. dx(3) ) then
     print *, "umax = ", umax, ", vmax = ", vmax, ", wmax = ", wmax, ", dt = ", dt, ", dx = ", dx
     call bl_error("CFL violation. Use smaller adv.cfl.")
  end if

  ! call a function to compute flux
  call mol2ndord_flux_3d(lo, hi, dt, dx, &
                         uin, ui_lo, ui_hi, &
                         vx, vx_lo, vx_hi, &
                         vy, vy_lo, vy_hi, &
                         vz, vz_lo, vz_hi, &
                         flxx, fx_lo, fx_hi, &
                         flxy, fy_lo, fy_hi, &
                         flxz, fz_lo, fz_hi, &
                         phix, phiy, phiz,  &
                         slope, glo, ghi,nu, uselimit)

  ! Do a conservative update
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           !notice that some mad scientist has reversed the flux difference
           ! rather than just put a negative sign out front.
           dphdtout(i,j,k) =  &
                ( (flxx(i,j,k) - flxx(i+1,j  ,k  )) /dx(1) &
                + (flxy(i,j,k) - flxy(i  ,j+1,k  )) /dx(2) &
                + (flxz(i,j,k) - flxz(i  ,j,  k+1)) /dx(3) )
        enddo
     enddo
  enddo

  ! Scale by face area in order to correctly reflux because flux register requires this
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
           flxx(i,j,k) = flxx(i,j,k) * (dt * dx(2)*dx(3))
        enddo
     enddo
  enddo
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)+1 
        do i = lo(1), hi(1)
           flxy(i,j,k) = flxy(i,j,k) * (dt * dx(1)*dx(3))
        enddo
     enddo
  enddo
  do       k = lo(3), hi(3)+1
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           flxz(i,j,k) = flxz(i,j,k) * (dt * dx(1)*dx(2))
        enddo
     enddo
  enddo

  call bl_deallocate(phix  )
  call bl_deallocate(phiy  )
  call bl_deallocate(phiz  )
  call bl_deallocate(slope)

end subroutine advectDiffMOL2ndOrd


subroutine advectDiffMOL4thOrd(time, lo, hi, &
                 uin , ui_lo, ui_hi, &
                 dphidtout, uo_lo, uo_hi, &
                 vx  , vx_lo, vx_hi, &
                 vy  , vy_lo, vy_hi, &
                 vz  , vz_lo, vz_hi, &
                 flxx, fx_lo, fx_hi, &
                 flxy, fy_lo, fy_hi, &
                 flxz, fz_lo, fz_hi, &
                 dx,dt,nu, &
                 deblocell, debhicell, &
                 hisidedebfacelo, hisidedebfacehi, &
                 losidedebfacelo, losidedebfacehi, printstuff, uselimit &
     ) bind(C, name="advectDiffMOL4thOrd")
  
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use compute_flux_module, only : mol4thord_flux_3d_nolimit, mol4thord_flux_3d_limited
  implicit none

  integer, intent(in) :: lo(3), hi(3), printstuff, uselimit
  integer, intent(in) :: deblocell(3), debhicell(3)
  integer, intent(in) :: hisidedebfacelo(3), hisidedebfacehi(3)
  integer, intent(in) :: losidedebfacelo(3), losidedebfacehi(3)
  double precision, intent(in) :: dx(3), dt, time, nu
  integer, intent(in) :: ui_lo(3), ui_hi(3)
  integer, intent(in) :: uo_lo(3), uo_hi(3)
  integer, intent(in) :: vx_lo(3), vx_hi(3)
  integer, intent(in) :: vy_lo(3), vy_hi(3)
  integer, intent(in) :: vz_lo(3), vz_hi(3)
  integer, intent(in) :: fx_lo(3), fx_hi(3)
  integer, intent(in) :: fy_lo(3), fy_hi(3)
  integer, intent(in) :: fz_lo(3), fz_hi(3)
  double precision, intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3))
  double precision, intent(inout) :: dphidtout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3))
  double precision, intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2),vx_lo(3):vx_hi(3))
  double precision, intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2),vy_lo(3):vy_hi(3))
  double precision, intent(in   ) :: vz  (vz_lo(1):vz_hi(1),vz_lo(2):vz_hi(2),vz_lo(3):vz_hi(3))
  double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
  double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
  double precision, intent(  out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))

  integer :: i, j, k
  integer :: glo(3), ghi(3)
  double precision ::  umax, vmax, wmax
! integer ::  numphi
! double precision::  phitot

  ! Some compiler may not support 'contiguous'.  Remove it in that case.
  double precision, dimension(:,:,:), pointer, contiguous :: &
         fluxptx, fluxpty, fluxptz,  phiptx, phipty, phiptz, phiavex, phiavey, phiavez, phiptcc



  glo = lo - 3
  ghi = hi + 3

  call bl_allocate(fluxptx  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(fluxpty  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(fluxptz  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate( phiptx  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate( phipty  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate( phiptz  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiavex  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiavey  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiavez  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))
  call bl_allocate(phiptcc  ,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))

  
  ! We like to allocate these **pointers** here and then pass them to a function
  ! to remove their pointerness for performance, because normally pointers could
  ! be aliasing.  We need to use pointers instead of allocatable arrays because
  ! we like to use BoxLib's bl_allocate to allocate memeory instead of the intrinsic
  ! allocate.  Bl_allocate is much faster than allocate inside OMP.  
  ! Note that one MUST CALL BL_DEALLOCATE.

  ! check if CFL condition is violated.
  umax = maxval(abs(vx))
  vmax = maxval(abs(vy))
  wmax = maxval(abs(vz))
  if ( umax*dt .ge. dx(1) .or. &
       vmax*dt .ge. dx(2) .or. &
       wmax*dt .ge. dx(3) ) then
     print *, "umax = ", umax, ", vmax = ", vmax, ", wmax = ", wmax, ", dt = ", dt, ", dx = ", dx
     call bl_error("CFL violation. Use smaller adv.cfl.")
  end if

  ! call a function to compute flux
  if(uselimit.eq.1) then
     call mol4thord_flux_3d_limited(&
          lo, hi, dt, dx, &
          uin, ui_lo, ui_hi, &
          vx, vx_lo, vx_hi, &
          vy, vy_lo, vy_hi, &
          vz, vz_lo, vz_hi, &
          flxx, fx_lo, fx_hi, &
          flxy, fy_lo, fy_hi, &
          flxz, fz_lo, fz_hi, &
          fluxptx, phiptx, phiavex,&
          fluxpty, phipty, phiavey,&
          fluxptz, phiptz, phiavez,&
          phiptcc, glo, ghi,nu, &
          deblocell, debhicell, &
          hisidedebfacelo, hisidedebfacehi,&
          losidedebfacelo, losidedebfacehi, printstuff)
  else
     call mol4thord_flux_3d_nolimit(&
          lo, hi, dt, dx, &
          uin, ui_lo, ui_hi, &
          vx, vx_lo, vx_hi, &
          vy, vy_lo, vy_hi, &
          vz, vz_lo, vz_hi, &
          flxx, fx_lo, fx_hi, &
          flxy, fy_lo, fy_hi, &
          flxz, fz_lo, fz_hi, &
          fluxptx, phiptx, phiavex,&
          fluxpty, phipty, phiavey,&
          fluxptz, phiptz, phiavez,&
          phiptcc, glo, ghi,nu, &
          deblocell, debhicell, &
          hisidedebfacelo, hisidedebfacehi,&
          losidedebfacelo, losidedebfacehi, printstuff)
  endif
  ! Do a conservative update
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           !notice that some mad scientist has reversed the flux difference
           ! rather than just put a negative sign out front.
           dphidtout(i,j,k) =  &
                ( (flxx(i,j,k) - flxx(i+1,j  ,k  )) /dx(1) &
                + (flxy(i,j,k) - flxy(i  ,j+1,k  )) /dx(2) &
                + (flxz(i,j,k) - flxz(i  ,j,  k+1)) /dx(3) )
        enddo
     enddo
  enddo

!  numphi = 0
!  phitot = 0.0d0
!  do       k = hisidedebfacelo(3), hisidedebfacehi(3)
!     do    j = hisidedebfacelo(2), hisidedebfacehi(2)
!        do i = hisidedebfacelo(1), hisidedebfacehi(1)
!           numphi = numphi + 1
!           phitot = phitot +flxx(i,j,k)
!           !          print*, "*** i j phiave = ", i, j, phiavex(i,j), "****"
!        enddo
!     enddo
!  enddo
!  
!  numphi = 0
!  phitot = 0.0d0
!  do       k = losidedebfacelo(3), losidedebfacehi(3)
!     do    j = losidedebfacelo(2), losidedebfacehi(2)
!        do i = losidedebfacelo(1), losidedebfacehi(1)
!           numphi = numphi + 1
!           phitot = phitot +flxx(i,j,k)
!           !          print*, "*** i j phiave = ", i, j, phiavex(i,j), "****"
!        enddo
!     enddo
!  enddo
!  if(numphi .gt. 0) then
!     print*, "**************** numphi, final lo x flux = ", numphi, phitot/numphi
!  endif
!  
!  numphi = 0
!  phitot = 0.0d0
!  do       k = deblocell(3), debhicell(3)
!     do    j = deblocell(2), debhicell(2)
!        do i = deblocell(1), debhicell(1)
!           numphi = numphi + 1
!           phitot = phitot + dphidtout(i,j,k)
!           !          print*, "*** i j phiave = ", i, j, phiavex(i,j), "****"
!        enddo
!     enddo
!  enddo
!  
!  if(numphi .gt. 0) then
!     print*, "**************** ndphidt, final dphidt = ", numphi, phitot/numphi
!  endif

  ! Scale by face area in order to correctly reflux because flux register requires this
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
           flxx(i,j,k) = flxx(i,j,k) * (dt * dx(2)*dx(3))
        enddo
     enddo
  enddo
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)+1 
        do i = lo(1), hi(1)
           flxy(i,j,k) = flxy(i,j,k) * (dt * dx(1)*dx(3))
        enddo
     enddo
  enddo
  do       k = lo(3), hi(3)+1
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           flxz(i,j,k) = flxz(i,j,k) * (dt * dx(1)*dx(2))
        enddo
     enddo
  enddo

  call bl_deallocate(fluxptx)
  call bl_deallocate(fluxpty)
  call bl_deallocate(fluxptz)
  call bl_deallocate( phiptx)
  call bl_deallocate( phipty)
  call bl_deallocate( phiptz)
  call bl_deallocate(phiavex)
  call bl_deallocate(phiavey)
  call bl_deallocate(phiavez)
  call bl_deallocate(phiptcc)


end subroutine advectDiffMOL4thOrd
