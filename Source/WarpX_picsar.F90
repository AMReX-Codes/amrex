#if (BL_SPACEDIM == 3)

#define WRPX_PXR_GETEB_ENERGY_CONSERVING geteb3d_energy_conserving_generic
#define WRPX_PXR_CURRENT_DEPOSITION depose_jxjyjz_generic

#elif (BL_SPACEDIM == 2)

#define WRPX_PXR_GETEB_ENERGY_CONSERVING geteb2dxz_energy_conserving
#define WRPX_PXR_CURRENT_DEPOSITION depose_jxjyjz_generic_2d

#endif

#define LVECT_CURRDEPO 8_c_long
#define LVECT_FIELDGATHE 64_c_long

! _________________________________________________________________
!
!> @brief
!> Module that contains subroutines to be called with AMReX
!> and that uses subroutines of Picsar
!>
!> @details
!> This avoids the use of interface with bind in the core of Picsar
!> This enables the use of integer in AMReX and Logical in Picsar
!> wihtout compatibility issue
!>
!> @author
!> Weiqun Zhang
!> Ann Almgren
!> Remi Lehe
!> Mathieu Lobet
!>
module warpx_to_pxr_module
! _________________________________________________________________

  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  use constants

  implicit none

  integer, parameter :: pxr_logical = 8

contains

  ! _________________________________________________________________
  !>
  !> @brief
  !> Main subroutine for the field gathering process
  !>
  !> @param[in] np number of particles
  !> @param[in] xp,yp,zp particle position arrays
  !> @param[in] ex,ey,ez particle electric fields in each direction
  !> @param[in] bx,by,bz particle magnetic fields in each direction
  !> @param[in] xmin,ymin,zmin tile grid minimum position
  !> @param[in] dx,dy,dz space discretization steps
  !> @param[in] nox,noy,noz interpolation order
  !> @param[in] exg,eyg,ezg electric field grid arrays
  !> @param[in] bxg,byg,bzg electric field grid arrays
  !> @param[in] lvect vector length
  !>
  subroutine warpx_geteb_energy_conserving(np,xp,yp,zp, &
       ex,ey,ez,bx,by,bz,xmin,ymin,zmin,dx,dy,dz,nox,noy,noz, &
       exg,exg_ng,exg_ntot,eyg,eyg_ng,eyg_ntot,ezg,ezg_ng,ezg_ntot, &
       bxg,bxg_ng,bxg_ntot,byg,byg_ng,byg_ntot,bzg,bzg_ng,bzg_ntot, &
       ll4symtry,l_lower_order_in_v, &
       lvect,&
       field_gathe_algo) &
       bind(C, name="warpx_geteb_energy_conserving")

    integer, intent(in) :: exg_ntot(BL_SPACEDIM), eyg_ntot(BL_SPACEDIM), ezg_ntot(BL_SPACEDIM), &
                           bxg_ntot(BL_SPACEDIM), byg_ntot(BL_SPACEDIM), bzg_ntot(BL_SPACEDIM)
    integer(c_long), intent(in) :: exg_ng, eyg_ng, ezg_ng, bxg_ng, byg_ng, bzg_ng
    integer(c_long), intent(in) :: field_gathe_algo
    integer(c_long), intent(in) :: np,nox,noy,noz
    integer(c_int), intent(in)  :: ll4symtry,l_lower_order_in_v
    integer(c_long),intent(in)   :: lvect
    real(amrex_real), intent(in), dimension(np) :: xp,yp,zp
    real(amrex_real), intent(out), dimension(np) :: ex,ey,ez,bx,by,bz
    real(amrex_real), intent(in)    :: xmin,ymin,zmin,dx,dy,dz
    real(amrex_real),intent(in):: exg(*), eyg(*), ezg(*), bxg(*), byg(*), bzg(*)
    logical(pxr_logical) :: pxr_ll4symtry, pxr_l_lower_order_in_v

    ! Compute the number of valid cells and guard cells
    integer(c_long) :: exg_nvalid(3), eyg_nvalid(3), ezg_nvalid(3),    &
                       bxg_nvalid(3), byg_nvalid(3), bzg_nvalid(3),    &
                       exg_nguards(3), eyg_nguards(3), ezg_nguards(3), &
                       bxg_nguards(3), byg_nguards(3), bzg_nguards(3)
#if (BL_SPACEDIM == 3)
    exg_nvalid = exg_ntot - 2*exg_ng
    eyg_nvalid = eyg_ntot - 2*eyg_ng
    ezg_nvalid = ezg_ntot - 2*ezg_ng
    bxg_nvalid = bxg_ntot - 2*bxg_ng
    byg_nvalid = byg_ntot - 2*byg_ng
    bzg_nvalid = bzg_ntot - 2*bzg_ng
    exg_nguards = exg_ng
    eyg_nguards = eyg_ng
    ezg_nguards = ezg_ng
    bxg_nguards = bxg_ng
    byg_nguards = byg_ng
    bzg_nguards = bzg_ng
#elif (BL_SPACEDIM == 2)
    ! The first and last index of `exg_ntot` correspond to the first and last index of the `exg_nvalid`,
    ! even though they do not have the same length
    ! Set the middle index of `exg_nvalid` to 1
    exg_nvalid(1:3:2) = exg_ntot(1:2) - 2*exg_ng
    exg_nvalid(2) = 1
    eyg_nvalid(1:3:2) = eyg_ntot(1:2) - 2*eyg_ng
    eyg_nvalid(2) = 1
    ezg_nvalid(1:3:2) = ezg_ntot(1:2) - 2*ezg_ng
    ezg_nvalid(2) = 1
    bxg_nvalid(1:3:2) = bxg_ntot(1:2) - 2*bxg_ng
    bxg_nvalid(2) = 1
    byg_nvalid(1:3:2) = byg_ntot(1:2) - 2*byg_ng
    byg_nvalid(2) = 1
    bzg_nvalid(1:3:2) = bzg_ntot(1:2) - 2*bzg_ng
    bzg_nvalid(2) = 1
    ! Set the number of guard cells in the middle index to 0
    exg_nguards = (/ exg_ng, 0_c_long, exg_ng /)
    eyg_nguards = (/ eyg_ng, 0_c_long, eyg_ng /)
    ezg_nguards = (/ ezg_ng, 0_c_long, ezg_ng /)
    bxg_nguards = (/ bxg_ng, 0_c_long, bxg_ng /)
    byg_nguards = (/ byg_ng, 0_c_long, byg_ng /)
    bzg_nguards = (/ bzg_ng, 0_c_long, bzg_ng /)
#endif

    pxr_ll4symtry = ll4symtry .eq. 1
    pxr_l_lower_order_in_v = l_lower_order_in_v .eq. 1

    CALL WRPX_PXR_GETEB_ENERGY_CONSERVING(np,xp,yp,zp, &
         ex,ey,ez,bx,by,bz,xmin,ymin,zmin,dx,dy,dz,nox,noy,noz, &
         exg,exg_nguards,exg_nvalid,&
         eyg,eyg_nguards,eyg_nvalid,&
         ezg,ezg_nguards,ezg_nvalid,&
         bxg,bxg_nguards,bxg_nvalid,&
         byg,byg_nguards,byg_nvalid,&
         bzg,bzg_nguards,bzg_nvalid,&
	 pxr_ll4symtry, pxr_l_lower_order_in_v, &
	 lvect, field_gathe_algo )

  end subroutine warpx_geteb_energy_conserving

! _________________________________________________________________
!>
!> @brief
!> Main subroutine for the charge deposition
!>
!> @details
!> This subroutines enable to controle the interpolation order
!> via the parameters nox,noy,noz and the type of algorithm via
!> the parameter charge_depo_algo
!
!> @param[inout] rho charge array
!> @param[in] np number of particles
!> @param[in] xp,yp,zp particle position arrays
!> @param[in] w particle weight arrays
!> @param[in] q particle species charge
!> @param[in] xmin,ymin,zmin tile grid minimum position
!> @param[in] dx,dy,dz space discretization steps
!> @param[in] nx,ny,nz number of cells
!> @param[in] nxguard,nyguard,nzguard number of guard cells
!> @param[in] nox,noy,noz interpolation order
!> @param[in] lvect vector length
!> @param[in] charge_depo_algo algorithm choice for the charge deposition
!>
subroutine warpx_charge_deposition(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
   nxguard,nyguard,nzguard,nox,noy,noz,lvect,charge_depo_algo) &
  bind(C, name="warpx_charge_deposition")

  integer(c_long), intent(IN)                                  :: np
  integer(c_long), intent(IN)                                  :: nx,ny,nz
  integer(c_long), intent(IN)                                  :: nxguard,nyguard,nzguard
  integer(c_long), intent(IN)                                  :: nox,noy,noz
  real(amrex_real), intent(IN OUT), dimension(-nxguard:nx+nxguard,&
       &                                  -nyguard:ny+nyguard,&
       &                                  -nzguard:nz+nzguard) :: rho
  real(amrex_real), intent(IN)                                     :: q
  real(amrex_real), intent(IN)                                     :: dx,dy,dz
  real(amrex_real), intent(IN)                                     :: xmin,ymin,zmin
  real(amrex_real), intent(IN),  dimension(np)                     :: xp,yp,zp,w
  integer(c_long), intent(IN)                                  :: lvect
  integer(c_long), intent(IN)                                  :: charge_depo_algo


  ! Dimension 3
#if (BL_SPACEDIM==3)

  SELECT CASE(charge_depo_algo)

  ! Scalar classical charge deposition subroutines
  CASE(1)
    IF ((nox.eq.1).and.(noy.eq.1).and.(noz.eq.1)) THEN

      CALL depose_rho_scalar_1_1_1(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
    nxguard,nyguard,nzguard,lvect)

    ELSE IF ((nox.eq.2).and.(noy.eq.2).and.(noz.eq.2)) THEN

      CALL depose_rho_scalar_2_2_2(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
    nxguard,nyguard,nzguard,lvect)

    ELSE IF ((nox.eq.3).and.(noy.eq.3).and.(noz.eq.3)) THEN

      CALL depose_rho_scalar_3_3_3(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
    nxguard,nyguard,nzguard,lvect)

    ELSE
      CALL pxr_depose_rho_n(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
                  nxguard,nyguard,nzguard,nox,noy,noz, &
                  .TRUE._c_long,.FALSE._c_long)
    ENDIF

  ! Optimized subroutines
  CASE DEFAULT

    IF ((nox.eq.1).and.(noy.eq.1).and.(noz.eq.1)) THEN
      CALL depose_rho_vecHVv2_1_1_1(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
               nxguard,nyguard,nzguard,lvect)
    ELSE
      CALL pxr_depose_rho_n(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
                  nxguard,nyguard,nzguard,nox,noy,noz, &
                  .TRUE._c_long,.FALSE._c_long)
    ENDIF
  END SELECT

  ! Dimension 2
#elif (BL_SPACEDIM==2)


#endif

 end subroutine warpx_charge_deposition

  ! _________________________________________________________________
  !>
  !> @brief
  !> Main subroutine for the current deposition
  !>
  !> @details
  !> This subroutines enable to controle the interpolation order
  !> via the parameters nox,noy,noz and the type of algorithm via
  !> the parameter current_depo_algo
  !
  !> @param[inout] jx,jy,jz current arrays
  !> @param[in] np number of particles
  !> @param[in] xp,yp,zp particle position arrays
  !> @param[in] uxp,uyp,uzp particle momentum arrays
  !> @param[in] gaminv inverve of the particle Lorentz factor (array)
  !> @param[in] w particle weight arrays
  !> @param[in] q particle species charge
  !> @param[in] xmin,ymin,zmin tile grid minimum position
  !> @param[in] dx,dy,dz space discretization steps
  !> @param[in] nx,ny,nz number of cells
  !> @param[in] nxguard,nyguard,nzguard number of guard cells
  !> @param[in] nox,noy,noz interpolation order
  !> @param[in] lvect vector length
  !> @param[in] charge_depo_algo algorithm choice for the charge deposition
  !>
  subroutine warpx_current_deposition( &
    jx,jx_ng,jx_ntot,jy,jy_ng,jy_ntot,jz,jz_ng,jz_ntot, &
    np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,ymin,zmin,dt,dx,dy,dz,nox,noy,noz,&
    lvect,current_depo_algo) &
    bind(C, name="warpx_current_deposition")

    integer, intent(in) :: jx_ntot(3), jy_ntot(3), jz_ntot(3)
    integer(c_long), intent(in) :: jx_ng, jy_ng, jz_ng
    integer(c_long), intent(IN)                                  :: np
    integer(c_long), intent(IN)                                  :: nox,noy,noz

    real(amrex_real), intent(IN OUT):: jx(*), jy(*), jz(*)
    real(amrex_real), intent(IN)                                     :: q
    real(amrex_real), intent(IN)                                     :: dx,dy,dz
    real(amrex_real), intent(IN)                                     :: dt
    real(amrex_real), intent(IN)                                     :: xmin,ymin,zmin
    real(amrex_real), intent(IN),  dimension(np)                     :: xp,yp,zp,w
    real(amrex_real), intent(IN),  dimension(np)                     :: uxp,uyp,uzp
    real(amrex_real), intent(IN),  dimension(np)                     :: gaminv
    integer(c_long), intent(IN)                                   :: lvect
    integer(c_long), intent(IN)                                   :: current_depo_algo

    ! Maintain variables nx, ny, nz, nxguard, nyguard, nzguard for compilation
    ! although they will not be used in the future
    integer(c_long) :: nx=0, ny=0, nz=0
    integer(c_long) :: nxguard=0, nyguard=0, nzguard=0

    ! Compute the number of valid cells and guard cells
    integer(c_long) :: jx_nvalid(3), jy_nvalid(3), jz_nvalid(3), &
                       jx_nguards(3), jy_nguards(3), jz_nguards(3)
#if (BL_SPACEDIM == 3)
   jx_nvalid = jx_ntot - 2*jx_ng
   jy_nvalid = jy_ntot - 2*jy_ng
   jz_nvalid = jz_ntot - 2*jz_ng
   jx_nguards = jx_ng
   jy_nguards = jy_ng
   jz_nguards = jz_ng
#elif (BL_SPACEDIM == 2)
   ! The first and last index of `jx_ntot` correspond to the first and last index of the `jx_nvalid`,
   ! even though they do not have the same length
   ! Set the middle index of `jx_nvalid` to 1
   jx_nvalid(1:3:2) = jx_ntot(1:2) - 2*jx_ng
   jx_nvalid(2) = 1
   jy_nvalid(1:3:2) = jy_ntot(1:2) - 2*jy_ng
   jy_nvalid(2) = 1
   jz_nvalid(1:3:2) = jz_ntot(1:2) - 2*jz_ng
   jz_nvalid(2) = 1
   ! Set the number of guard cells in the middle index to 0
   jx_nguards = (/ jx_ng, 0_c_long, jx_ng /)
   jy_nguards = (/ jy_ng, 0_c_long, jy_ng /)
   jz_nguards = (/ jz_ng, 0_c_long, jz_ng /)
#endif

! Dimension 3
#if (BL_SPACEDIM==3)
   CALL WRPX_PXR_CURRENT_DEPOSITION(        &
        jx,jx_nguards,jx_nvalid,            &
        jy,jy_nguards,jy_nvalid,            &
        jz,jz_nguards,jz_nvalid,            &
        np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q, &
        xmin,ymin,zmin,dt,dx,dy,dz,         &
        nox,noy,noz,current_depo_algo)
! Dimension 2
#elif (BL_SPACEDIM==2)
        CALL WRPX_PXR_CURRENT_DEPOSITION(jx,jy,jz,np, &
        xp,yp,zp,uxp,uyp,uzp,gaminv,w,q,xmin,zmin, &
        dt,dx,dz,nx,nz,nxguard,nzguard,nox,noz,lvect)
#endif

  end subroutine

  ! _________________________________________________________________
  !>
  !> @brief
  !> Main subroutine for the particle pusher (velocity and position)
  !>
  !> @param[in] np number of super-particles
  !> @param[in] xp,yp,zp particle position arrays
  !> @param[in] uxp,uyp,uzp normalized momentum in each direction
  !> @param[in] gaminv particle Lorentz factors
  !> @param[in] ex,ey,ez particle electric fields in each direction
  !> @param[in] bx,by,bz particle magnetic fields in each direction
  !> @param[in] q charge
  !> @param[in] m masse
  !> @param[in] dt time step
  !> @param[in] particle_pusher_algo Particle pusher algorithm
  subroutine warpx_particle_pusher(np,xp,yp,zp,uxp,uyp,uzp, &
                                  gaminv,&
                                  ex,ey,ez,bx,by,bz,q,m,dt, &
                                  particle_pusher_algo) &
       bind(C, name="warpx_particle_pusher")

    INTEGER(c_long), INTENT(IN)   :: np
    REAL(amrex_real),INTENT(INOUT)    :: gaminv(np)
    REAL(amrex_real),INTENT(INOUT)    :: xp(np),yp(np),zp(np)
    REAL(amrex_real),INTENT(INOUT)    :: uxp(np),uyp(np),uzp(np)
    REAL(amrex_real),INTENT(IN)       :: ex(np),ey(np),ez(np)
    REAL(amrex_real),INTENT(IN)       :: bx(np),by(np),bz(np)
    REAL(amrex_real),INTENT(IN)       :: q,m,dt
    INTEGER(c_long), INTENT(IN)   :: particle_pusher_algo

    SELECT CASE (particle_pusher_algo)

    !! Vay pusher -- Full push
    CASE (1_c_long)
!#if (BL_SPACEDIM == 3)
      CALL pxr_ebcancelpush3d(np,uxp,uyp,uzp,gaminv, &
                                 ex,ey,ez,  &
                                 bx,by,bz,q,m,dt,0_c_long)
!#else
!      call bl_error("Is there a 2d Vay pusher implemented?")
!#endif
    CASE DEFAULT

      ! Momentum pusher in a single loop
      CALL pxr_boris_push_u_3d(np,uxp,uyp,uzp,&
                                     gaminv, &
                                     ex,ey,ez, &
                                     bx,by,bz, &
                                     q,m,dt)

      ! Momentum pusher by block
!       CALL pxr_boris_push_u_3d_block(np,uxp,uyp,uzp,&
!                                      gaminv, &
!                                      ex,ey,ez, &
!                                      bx,by,bz, &
!                                      q,m,dt,lvect)

      !! --- Push velocity with E half step
!       CALL pxr_epush_v(np,uxp,uyp,uzp,      &
!                       ex,ey,ez,q,m,dt*0.5_amrex_real)
      !! --- Set gamma of particles
!       CALL pxr_set_gamma(np,uxp,uyp,uzp,gaminv)
      !! --- Push velocity with B
!       CALL pxr_bpush_v(np,uxp,uyp,uzp,gaminv,      &
!                       bx,by,bz,q,m,dt)
      !!! --- Push velocity with E half step
!       CALL pxr_epush_v(np,uxp,uyp,uzp,      &
!                       ex,ey,ez,q,m,dt*0.5_amrex_real)
      !! --- Set gamma of particles
!       CALL pxr_set_gamma(np,uxp,uyp,uzp,gaminv)

    END SELECT

    !!!! --- push particle species positions a time step
#if (BL_SPACEDIM == 3)
    CALL pxr_pushxyz(np,xp,yp,zp,uxp,uyp,uzp,gaminv,dt)
#elif (BL_SPACEDIM == 2)
    CALL pxr_pushxz(np,xp,zp,uxp,uzp,gaminv,dt)
#endif

  end subroutine

  ! _________________________________________________________________
  !>
  !> @brief
  !> Main subroutine for the particle pusher of positions
  !>
  !> @param[in] np number of super-particles
  !> @param[in] xp,yp,zp particle position arrays
  !> @param[in] uxp,uyp,uzp normalized momentum in each direction
  !> @param[in] gaminv particle Lorentz factors
  !> @param[in] dt time step
  !> @param[in] particle_pusher_algo Particle pusher algorithm
  subroutine warpx_particle_pusher_positions(np,xp,yp,zp,uxp,uyp,uzp, &
                                  gaminv,dt) &
       bind(C, name="warpx_particle_pusher_positions")

    INTEGER(c_long), INTENT(IN)   :: np
    REAL(amrex_real),INTENT(INOUT)    :: gaminv(np)
    REAL(amrex_real),INTENT(INOUT)    :: xp(np),yp(np),zp(np)
    REAL(amrex_real),INTENT(INOUT)    :: uxp(np),uyp(np),uzp(np)
    REAL(amrex_real),INTENT(IN)       :: dt

    CALL pxr_set_gamma(np,uxp,uyp,uzp,gaminv)

    !!!! --- push particle species positions a time step
#if (BL_SPACEDIM == 3)
    CALL pxr_pushxyz(np,xp,yp,zp,uxp,uyp,uzp,gaminv,dt)
#elif (BL_SPACEDIM == 2)
    CALL pxr_pushxz(np,xp,zp,uxp,uzp,gaminv,dt)
#endif

  end subroutine

end module warpx_to_pxr_module
