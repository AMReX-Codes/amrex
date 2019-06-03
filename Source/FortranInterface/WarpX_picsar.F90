#if (AMREX_SPACEDIM == 3)

#define WRPX_PXR_GETEB_ENERGY_CONSERVING  geteb3d_energy_conserving_generic
#define WRPX_PXR_CURRENT_DEPOSITION       depose_jxjyjz_generic
#define WRPX_PXR_PUSH_BVEC                pxrpush_em3d_bvec
#define WRPX_PXR_PUSH_BVEC_CKC            pxrpush_em3d_bvec_ckc
#define WRPX_PXR_PUSH_EVEC_F              pxrpush_em3d_evec_f
#define WRPX_PXR_PUSH_EVEC_F_CKC          pxrpush_em3d_evec_f_ckc
#define WRPX_PXR_PUSH_EVEC                pxrpush_em3d_evec

#elif (AMREX_SPACEDIM == 2)

#define WRPX_PXR_PUSH_BVEC_CKC            pxrpush_em2d_bvec_ckc
#define WRPX_PXR_PUSH_EVEC_F              pxrpush_em2d_evec_f
#define WRPX_PXR_PUSH_EVEC_F_CKC          pxrpush_em2d_evec_f_ckc

#ifdef WARPX_RZ

#define WRPX_PXR_GETEB_ENERGY_CONSERVING  geteb2drz_energy_conserving_generic
#define WRPX_PXR_CURRENT_DEPOSITION       depose_jrjtjz_generic_rz
#define WRPX_PXR_RZ_VOLUME_SCALING_RHO    apply_rz_volume_scaling_rho
#define WRPX_PXR_RZ_VOLUME_SCALING_J      apply_rz_volume_scaling_j
#define WRPX_PXR_PUSH_BVEC                pxrpush_emrz_bvec
#define WRPX_PXR_PUSH_EVEC                pxrpush_emrz_evec

#else

#define WRPX_PXR_GETEB_ENERGY_CONSERVING  geteb2dxz_energy_conserving_generic
#define WRPX_PXR_CURRENT_DEPOSITION       depose_jxjyjz_generic_2d
#define WRPX_PXR_PUSH_BVEC                pxrpush_em2d_bvec
#define WRPX_PXR_PUSH_EVEC                pxrpush_em2d_evec

#endif

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
  !> @param[in] ixyzmin tile grid minimum index
  !> @param[in] xmin,ymin,zmin tile grid minimum position
  !> @param[in] dx,dy,dz space discretization steps
  !> @param[in] xyzmin grid minimum position
  !> @param[in] dxyz space discretization steps
  !> @param[in] nox,noy,noz interpolation order
  !> @param[in] exg,eyg,ezg electric field grid arrays
  !> @param[in] bxg,byg,bzg electric field grid arrays
  !> @param[in] lvect vector length
  !>
  subroutine warpx_geteb_energy_conserving(np,xp,yp,zp, &
       ex,ey,ez,bx,by,bz,ixyzmin,xmin,ymin,zmin,dx,dy,dz,nox,noy,noz, &
       exg,exg_lo,exg_hi,eyg,eyg_lo,eyg_hi,ezg,ezg_lo,ezg_hi, &
       bxg,bxg_lo,bxg_hi,byg,byg_lo,byg_hi,bzg,bzg_lo,bzg_hi, &
       ll4symtry,l_lower_order_in_v, l_nodal,&
       lvect,field_gathe_algo) &
       bind(C, name="warpx_geteb_energy_conserving")

    integer, intent(in) :: exg_lo(AMREX_SPACEDIM), eyg_lo(AMREX_SPACEDIM), ezg_lo(AMREX_SPACEDIM), &
                           bxg_lo(AMREX_SPACEDIM), byg_lo(AMREX_SPACEDIM), bzg_lo(AMREX_SPACEDIM)
    integer, intent(in) :: exg_hi(AMREX_SPACEDIM), eyg_hi(AMREX_SPACEDIM), ezg_hi(AMREX_SPACEDIM), &
                           bxg_hi(AMREX_SPACEDIM), byg_hi(AMREX_SPACEDIM), bzg_hi(AMREX_SPACEDIM)
    integer, intent(in) :: ixyzmin(AMREX_SPACEDIM)
    real(amrex_real), intent(in) :: xmin,ymin,zmin,dx,dy,dz
    integer(c_long), intent(in) :: field_gathe_algo
    integer(c_long), intent(in) :: np,nox,noy,noz
    integer(c_int), intent(in)  :: ll4symtry,l_lower_order_in_v, l_nodal
    integer(c_long),intent(in)   :: lvect
    real(amrex_real), intent(in), dimension(np) :: xp,yp,zp
    real(amrex_real), intent(out), dimension(np) :: ex,ey,ez,bx,by,bz
    real(amrex_real),intent(in):: exg(*), eyg(*), ezg(*), bxg(*), byg(*), bzg(*)
    logical(pxr_logical) :: pxr_ll4symtry, pxr_l_lower_order_in_v, pxr_l_nodal

    ! Compute the number of valid cells and guard cells
    integer(c_long) :: exg_nvalid(AMREX_SPACEDIM), eyg_nvalid(AMREX_SPACEDIM), ezg_nvalid(AMREX_SPACEDIM),    &
                       bxg_nvalid(AMREX_SPACEDIM), byg_nvalid(AMREX_SPACEDIM), bzg_nvalid(AMREX_SPACEDIM),    &
                       exg_nguards(AMREX_SPACEDIM), eyg_nguards(AMREX_SPACEDIM), ezg_nguards(AMREX_SPACEDIM), &
                       bxg_nguards(AMREX_SPACEDIM), byg_nguards(AMREX_SPACEDIM), bzg_nguards(AMREX_SPACEDIM)

    pxr_ll4symtry = ll4symtry .eq. 1
    pxr_l_lower_order_in_v = l_lower_order_in_v .eq. 1
    pxr_l_nodal = l_nodal .eq. 1

    exg_nguards = ixyzmin - exg_lo
    eyg_nguards = ixyzmin - eyg_lo
    ezg_nguards = ixyzmin - ezg_lo
    bxg_nguards = ixyzmin - bxg_lo
    byg_nguards = ixyzmin - byg_lo
    bzg_nguards = ixyzmin - bzg_lo
    exg_nvalid = exg_lo + exg_hi - 2_c_long*ixyzmin + 1_c_long
    eyg_nvalid = eyg_lo + eyg_hi - 2_c_long*ixyzmin + 1_c_long
    ezg_nvalid = ezg_lo + ezg_hi - 2_c_long*ixyzmin + 1_c_long
    bxg_nvalid = bxg_lo + bxg_hi - 2_c_long*ixyzmin + 1_c_long
    byg_nvalid = byg_lo + byg_hi - 2_c_long*ixyzmin + 1_c_long
    bzg_nvalid = bzg_lo + bzg_hi - 2_c_long*ixyzmin + 1_c_long

    CALL WRPX_PXR_GETEB_ENERGY_CONSERVING(np,xp,yp,zp, &
         ex,ey,ez,bx,by,bz,xmin,ymin,zmin,dx,dy,dz,nox,noy,noz, &
         exg,exg_nguards,exg_nvalid,&
         eyg,eyg_nguards,eyg_nvalid,&
         ezg,ezg_nguards,ezg_nvalid,&
         bxg,bxg_nguards,bxg_nvalid,&
         byg,byg_nguards,byg_nvalid,&
         bzg,bzg_nguards,bzg_nvalid,&
	 pxr_ll4symtry, pxr_l_lower_order_in_v, pxr_l_nodal, &
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

  integer(c_long), intent(IN)                   :: np
  integer(c_long), intent(IN)                   :: nx,ny,nz
  integer(c_long), intent(IN)                   :: nxguard,nyguard,nzguard
  integer(c_long), intent(IN)                   :: nox,noy,noz
  real(amrex_real), intent(IN OUT)              :: rho(*)
  real(amrex_real), intent(IN)                  :: q
  real(amrex_real), intent(IN)                  :: dx,dy,dz
  real(amrex_real), intent(IN)                  :: xmin,ymin,zmin
  real(amrex_real), intent(IN),  dimension(np)  :: xp,yp,zp,w
  integer(c_long), intent(IN)                   :: lvect
  integer(c_long), intent(IN)                   :: charge_depo_algo


  ! Dimension 3
#if (AMREX_SPACEDIM==3)

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

    ELSE IF ((nox.eq.2).and.(noy.eq.2).and.(noz.eq.2)) THEN
      CALL depose_rho_vecHVv2_2_2_2(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
                 nxguard,nyguard,nzguard,lvect)

    ELSE
      CALL pxr_depose_rho_n(rho,np,xp,yp,zp,w,q,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,&
                  nxguard,nyguard,nzguard,nox,noy,noz, &
                  .TRUE._c_long,.FALSE._c_long)
    ENDIF
  END SELECT

  ! Dimension 2
#elif (AMREX_SPACEDIM==2)

#ifdef WARPX_RZ
  logical(pxr_logical) :: l_2drz = .TRUE._c_long
#else
  logical(pxr_logical) :: l_2drz = .FALSE._c_long
#endif

  CALL pxr_depose_rho_n_2dxz(rho,np,xp,yp,zp,w,q,xmin,zmin,dx,dz,nx,nz,&
       nxguard,nzguard,nox,noz, &
       .TRUE._c_long, .FALSE._c_long, l_2drz, 0_c_long)

#endif

 end subroutine warpx_charge_deposition

  ! _________________________________________________________________
  !>
  !> @brief
  !> Applies the inverse volume scaling for RZ charge deposition
  !>
  !> @details
  !> The scaling is done for both single mode (FDTD) and
  !> multi mode (spectral) (todo)
  !
  !> @param[inout] rho charge array
  !> @param[in] rmin tile grid minimum radius
  !> @param[in] dr radial space discretization steps
  !> @param[in] nx,ny,nz number of cells
  !> @param[in] nxguard,nyguard,nzguard number of guard cells
  !>
  subroutine warpx_charge_deposition_rz_volume_scaling(rho,rho_ng,rho_ntot,rmin,dr) &
    bind(C, name="warpx_charge_deposition_rz_volume_scaling")

    integer, intent(in) :: rho_ntot(AMREX_SPACEDIM)
    integer(c_long), intent(in) :: rho_ng
    real(amrex_real), intent(IN OUT):: rho(*)
    real(amrex_real), intent(IN) :: rmin, dr

#ifdef WARPX_RZ
    integer(c_long) :: type_rz_depose = 1
#endif

    ! Compute the number of valid cells and guard cells
    integer(c_long) :: rho_nvalid(AMREX_SPACEDIM), rho_nguards(AMREX_SPACEDIM)
    rho_nvalid = rho_ntot - 2*rho_ng
    rho_nguards = rho_ng

#ifdef WARPX_RZ
    CALL WRPX_PXR_RZ_VOLUME_SCALING_RHO(   &
                 rho,rho_nguards,rho_nvalid, &
                 rmin,dr,type_rz_depose)
#endif

  end subroutine warpx_charge_deposition_rz_volume_scaling

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

    integer, intent(in) :: jx_ntot(AMREX_SPACEDIM), jy_ntot(AMREX_SPACEDIM), jz_ntot(AMREX_SPACEDIM)
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

    ! Compute the number of valid cells and guard cells
    integer(c_long) :: jx_nvalid(AMREX_SPACEDIM), jy_nvalid(AMREX_SPACEDIM), jz_nvalid(AMREX_SPACEDIM), &
                       jx_nguards(AMREX_SPACEDIM), jy_nguards(AMREX_SPACEDIM), jz_nguards(AMREX_SPACEDIM)
    jx_nvalid = jx_ntot - 2*jx_ng
    jy_nvalid = jy_ntot - 2*jy_ng
    jz_nvalid = jz_ntot - 2*jz_ng
    jx_nguards = jx_ng
    jy_nguards = jy_ng
    jz_nguards = jz_ng

! Dimension 3
#if (AMREX_SPACEDIM==3)
   CALL WRPX_PXR_CURRENT_DEPOSITION(        &
        jx,jx_nguards,jx_nvalid,            &
        jy,jy_nguards,jy_nvalid,            &
        jz,jz_nguards,jz_nvalid,            &
        np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q, &
        xmin,ymin,zmin,dt,dx,dy,dz,         &
        nox,noy,noz,current_depo_algo)
! Dimension 2
#elif (AMREX_SPACEDIM==2)
        CALL WRPX_PXR_CURRENT_DEPOSITION(   &
        jx,jx_nguards,jx_nvalid,            &
        jy,jy_nguards,jy_nvalid,            &
        jz,jz_nguards,jz_nvalid,            &
        np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q, &
        xmin,zmin,dt,dx,dz,nox,noz,lvect,   &
        current_depo_algo)
#endif

  end subroutine warpx_current_deposition

  ! _________________________________________________________________
  !>
  !> @brief
  !> Applies the inverse volume scaling for RZ current deposition
  !>
  !> @details
  !> The scaling is done for single mode only
  !
  !> @param[inout] jx,jy,jz current arrays
  !> @param[in] jx_ntot,jy_ntot,jz_ntot vectors with total number of
  !>            cells (including guard cells) along each axis for each current
  !> @param[in] jx_ng,jy_ng,jz_ng vectors with number of guard cells along each
  !>            axis for each current
  !> @param[in] rmin tile grid minimum radius
  !> @param[in] dr radial space discretization steps
  !>
  subroutine warpx_current_deposition_rz_volume_scaling( &
    jx,jx_ng,jx_ntot,jy,jy_ng,jy_ntot,jz,jz_ng,jz_ntot, &
    rmin,dr) &
    bind(C, name="warpx_current_deposition_rz_volume_scaling")

    integer, intent(in) :: jx_ntot(AMREX_SPACEDIM), jy_ntot(AMREX_SPACEDIM), jz_ntot(AMREX_SPACEDIM)
    integer(c_long), intent(in) :: jx_ng, jy_ng, jz_ng
    real(amrex_real), intent(IN OUT):: jx(*), jy(*), jz(*)
    real(amrex_real), intent(IN) :: rmin, dr

#ifdef WARPX_RZ
    integer(c_long) :: type_rz_depose = 1
#endif
    ! Compute the number of valid cells and guard cells
    integer(c_long) :: jx_nvalid(AMREX_SPACEDIM), jy_nvalid(AMREX_SPACEDIM), jz_nvalid(AMREX_SPACEDIM), &
                       jx_nguards(AMREX_SPACEDIM), jy_nguards(AMREX_SPACEDIM), jz_nguards(AMREX_SPACEDIM)
    jx_nvalid = jx_ntot - 2*jx_ng
    jy_nvalid = jy_ntot - 2*jy_ng
    jz_nvalid = jz_ntot - 2*jz_ng
    jx_nguards = jx_ng
    jy_nguards = jy_ng
    jz_nguards = jz_ng

#ifdef WARPX_RZ
    CALL WRPX_PXR_RZ_VOLUME_SCALING_J(   &
                 jx,jx_nguards,jx_nvalid, &
                 jy,jy_nguards,jy_nvalid, &
                 jz,jz_nguards,jz_nvalid, &
                 rmin,dr,type_rz_depose)
#endif

  end subroutine warpx_current_deposition_rz_volume_scaling

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
      CALL pxr_set_gamma(np,uxp,uyp,uzp,gaminv)

      CALL pxr_ebcancelpush3d(np,uxp,uyp,uzp,gaminv, &
                                 ex,ey,ez,  &
                                 bx,by,bz,q,m,dt,0_c_long)
    CASE DEFAULT

      ! Momentum pusher in a single loop
      CALL pxr_boris_push_u_3d(np,uxp,uyp,uzp,&
                                     gaminv, &
                                     ex,ey,ez, &
                                     bx,by,bz, &
                                     q,m,dt)

    END SELECT

    !!!! --- push particle species positions a time step
#if (AMREX_SPACEDIM == 3) || (defined WARPX_RZ)
    CALL pxr_pushxyz(np,xp,yp,zp,uxp,uyp,uzp,gaminv,dt)
#elif (AMREX_SPACEDIM == 2)
    CALL pxr_pushxz(np,xp,zp,uxp,uzp,gaminv,dt)
#endif

  end subroutine warpx_particle_pusher


  ! _________________________________________________________________
  !>
  !> @brief
  !> Main subroutine for the particle pusher (velocity)
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
  subroutine warpx_particle_pusher_momenta(np,xp,yp,zp,uxp,uyp,uzp, &
                                  gaminv,&
                                  ex,ey,ez,bx,by,bz,q,m,dt, &
                                  particle_pusher_algo) &
       bind(C, name="warpx_particle_pusher_momenta")

    INTEGER(c_long), INTENT(IN)   :: np
    REAL(amrex_real),INTENT(INOUT)    :: gaminv(np)
    REAL(amrex_real),INTENT(IN)       :: xp(np),yp(np),zp(np)
    REAL(amrex_real),INTENT(INOUT)    :: uxp(np),uyp(np),uzp(np)
    REAL(amrex_real),INTENT(IN)       :: ex(np),ey(np),ez(np)
    REAL(amrex_real),INTENT(IN)       :: bx(np),by(np),bz(np)
    REAL(amrex_real),INTENT(IN)       :: q,m,dt
    INTEGER(c_long), INTENT(IN)   :: particle_pusher_algo

    SELECT CASE (particle_pusher_algo)

    !! Vay pusher -- Full push
    CASE (1_c_long)
      CALL pxr_set_gamma(np,uxp,uyp,uzp,gaminv)

      CALL pxr_ebcancelpush3d(np,uxp,uyp,uzp,gaminv, &
                                 ex,ey,ez,  &
                                 bx,by,bz,q,m,dt,0_c_long)
    CASE DEFAULT

      ! Momentum pusher in a single loop
      CALL pxr_boris_push_u_3d(np,uxp,uyp,uzp,&
                                     gaminv, &
                                     ex,ey,ez, &
                                     bx,by,bz, &
                                     q,m,dt)

    END SELECT

  end subroutine warpx_particle_pusher_momenta

  ! _________________________________________________________________
  !>
  !> @brief
  !> Main subroutine for the Maxwell solver (E field)
  !>
  !> @param[in] xlo, xhi, ylo, yhi, zlo, zhi lower and higher bounds
  !>            where to apply the solver (PICSAR's valid cells)
  !> @param[inout] ex, ey, ez arrays for electric field to update
  !> @param[in] bx, by, bz, jx, jy, jz arrays for magnetic field and current
  !> @param[in] exlo, exhi, eylo, eyhi, ezlo, ezhi lower and higher bound for
  !>            the electric field arrays
  !> @param[in] bxlo, bxhi, bylo, byhi, bzlo, bzhi lower and higher bound for
  !>            the magnetic field arrays
  !> @param[in] jxlo, jxhi, jylo, jyhi, jzlo, jzhi lower and higher bound for
  !>            the current arrays
  !> @param[in] mudt normalized time step mu_0 * c**2 * dt
  !> @param[in] dtsdx, dtsdy, dtsdz factors c**2 * dt/(dx, dy, dz)
  subroutine warpx_push_evec( &
    xlo, xhi, ylo, yhi, zlo, zhi, &
    ex, exlo, exhi, &
    ey, eylo, eyhi, &
    ez, ezlo, ezhi, &
    bx, bxlo, bxhi, &
    by, bylo, byhi, &
    bz, bzlo, bzhi, &
    jx, jxlo, jxhi, &
    jy, jylo, jyhi, &
    jz, jzlo, jzhi, &
    mudt, dtsdx, dtsdy, dtsdz, xmin, dx) bind(C, name="warpx_push_evec")

    integer(c_int), intent(in) :: xlo(BL_SPACEDIM), xhi(BL_SPACEDIM), &
      ylo(BL_SPACEDIM), yhi(BL_SPACEDIM), zlo(BL_SPACEDIM), zhi(BL_SPACEDIM), &
      exlo(BL_SPACEDIM), exhi(BL_SPACEDIM), eylo(BL_SPACEDIM), eyhi(BL_SPACEDIM), &
      ezlo(BL_SPACEDIM), ezhi(BL_SPACEDIM), bxlo(BL_SPACEDIM), bxhi(BL_SPACEDIM), &
      bylo(BL_SPACEDIM), byhi(BL_SPACEDIM), bzlo(BL_SPACEDIM), bzhi(BL_SPACEDIM), &
      jxlo(BL_SPACEDIM), jxhi(BL_SPACEDIM), jylo(BL_SPACEDIM), jyhi(BL_SPACEDIM), &
      jzlo(BL_SPACEDIM), jzhi(BL_SPACEDIM)

    real(amrex_real), intent(IN OUT):: ex(*), ey(*), ez(*)

    real(amrex_real), intent(IN):: bx(*), by(*), bz(*), jx(*), jy(*), jz(*)

    real(amrex_real), intent(IN) :: mudt, dtsdx, dtsdy, dtsdz

    real(amrex_real), intent(IN) :: xmin, dx

    CALL WRPX_PXR_PUSH_EVEC(&
        xlo, xhi, ylo, yhi, zlo, zhi, &
        ex, exlo, exhi,&
        ey, eylo, eyhi,&
        ez, ezlo, ezhi,&
        bx, bxlo, bxhi,&
        by, bylo, byhi,&
        bz, bzlo, bzhi,&
        jx, jxlo, jxhi,&
        jy, jylo, jyhi,&
        jz, jzlo, jzhi,&
        mudt, dtsdx, dtsdy, dtsdz &
#ifdef WARPX_RZ
        ,xmin,dx &
#endif
        )

  end subroutine warpx_push_evec

  ! _________________________________________________________________
  !>
  !> @brief
  !> Main subroutine for the Maxwell solver (B field)
  !>
  !> @param[in] xlo, xhi, ylo, yhi, zlo, zhi lower and higher bounds
  !>            where to apply the solver (PICSAR's valid cells)
  !> @param[in] ex, ey, ez arrays for electric field
  !> @param[inout] bx, by, bz arrays for magnetic field to update
  !> @param[in] exlo, exhi, eylo, eyhi, ezlo, ezhi lower and higher bound for
  !>            the electric field arrays
  !> @param[in] bxlo, bxhi, bylo, byhi, bzlo, bzhi lower and higher bound for
  !>            the magnetic field arrays
  !> @param[in] dtsdx, dtsdy, dtsdz factors 0.5 * dt/(dx, dy, dz)
  subroutine warpx_push_bvec( &
    xlo, xhi, ylo, yhi, zlo, zhi, &
    ex, exlo, exhi, &
    ey, eylo, eyhi, &
    ez, ezlo, ezhi, &
    bx, bxlo, bxhi, &
    by, bylo, byhi, &
    bz, bzlo, bzhi, &
    dtsdx, dtsdy, dtsdz, xmin, dx, &
    maxwell_fdtd_solver_id) bind(C, name="warpx_push_bvec")

    integer(c_int), intent(in) :: xlo(BL_SPACEDIM), xhi(BL_SPACEDIM), &
      ylo(BL_SPACEDIM), yhi(BL_SPACEDIM), zlo(BL_SPACEDIM), zhi(BL_SPACEDIM), &
      exlo(BL_SPACEDIM), exhi(BL_SPACEDIM), eylo(BL_SPACEDIM), eyhi(BL_SPACEDIM), &
      ezlo(BL_SPACEDIM), ezhi(BL_SPACEDIM), bxlo(BL_SPACEDIM), bxhi(BL_SPACEDIM), &
      bylo(BL_SPACEDIM), byhi(BL_SPACEDIM), bzlo(BL_SPACEDIM), bzhi(BL_SPACEDIM), &
      maxwell_fdtd_solver_id

    real(amrex_real), intent(IN OUT):: ex(*), ey(*), ez(*)

    real(amrex_real), intent(IN):: bx(*), by(*), bz(*)

    real(amrex_real), intent(IN) :: dtsdx, dtsdy, dtsdz

    real(amrex_real), intent(IN) :: xmin, dx

    IF (maxwell_fdtd_solver_id .eq. 0) THEN
      ! Yee FDTD solver
      CALL WRPX_PXR_PUSH_BVEC( &
        xlo, xhi, ylo, yhi, zlo, zhi, &
      	ex, exlo, exhi, &
      	ey, eylo, eyhi, &
      	ez, ezlo, ezhi, &
      	bx, bxlo, bxhi, &
      	by, bylo, byhi, &
      	bz, bzlo, bzhi, &
      	dtsdx,dtsdy,dtsdz &
#ifdef WARPX_RZ
        ,xmin,dx &
#endif
        )
    ELSE IF (maxwell_fdtd_solver_id .eq. 1) THEN
      ! Cole-Karkkainen FDTD solver
      CALL WRPX_PXR_PUSH_BVEC_CKC( &
        xlo, xhi, ylo, yhi, zlo, zhi, &
      	ex, exlo, exhi, &
      	ey, eylo, eyhi, &
      	ez, ezlo, ezhi, &
      	bx, bxlo, bxhi, &
      	by, bylo, byhi, &
      	bz, bzlo, bzhi, &
      	dtsdx,dtsdy,dtsdz)
    ENDIF
  end subroutine warpx_push_bvec

  ! _________________________________________________________________
  !>
  !> @brief
  !> Subroutine for pushing E with the grad F terms
  !>
  !> @param[in] xlo, xhi, ylo, yhi, zlo, zhi lower and higher bounds
  !>            where to apply the solver (PICSAR's valid cells)
  !> @param[inout] ex, ey, ez arrays for electric field to update
  !> @param[in] f array for the charge-correction term
  !> @param[in] exlo, exhi, eylo, eyhi, ezlo, ezhi lower and higher bound for
  !>            the electric field arrays
  !> @param[in] flo, fhi, lower and higher bound for the charge-correction term
  !> @param[in] dtsdx_c2, dtsdy_c2, dtsdz_c2 factors dt/(dx, dy, dz)*c2
  subroutine warpx_push_evec_f( &
    xlo, xhi, ylo, yhi, zlo, zhi, &
    ex, exlo, exhi, &
    ey, eylo, eyhi, &
    ez, ezlo, ezhi, &
    f, flo, fhi, &
    dtsdx_c2, dtsdy_c2, dtsdz_c2, &
    maxwell_fdtd_solver_id) bind(C, name="warpx_push_evec_f")

    integer(c_int), intent(in) :: xlo(BL_SPACEDIM), xhi(BL_SPACEDIM), &
      ylo(BL_SPACEDIM), yhi(BL_SPACEDIM), zlo(BL_SPACEDIM), zhi(BL_SPACEDIM), &
      exlo(BL_SPACEDIM), exhi(BL_SPACEDIM), eylo(BL_SPACEDIM), eyhi(BL_SPACEDIM), &
      ezlo(BL_SPACEDIM), ezhi(BL_SPACEDIM), flo(BL_SPACEDIM), fhi(BL_SPACEDIM), &
      maxwell_fdtd_solver_id

    real(amrex_real), intent(IN OUT):: ex(*), ey(*), ez(*)

    real(amrex_real), intent(IN):: f

    real(amrex_real), intent(IN) :: dtsdx_c2, dtsdy_c2, dtsdz_c2

    IF (maxwell_fdtd_solver_id .eq. 0) THEN
      ! Yee FDTD solver
      CALL WRPX_PXR_PUSH_EVEC_F( &
        xlo, xhi, ylo, yhi, zlo, zhi, &
      	ex, exlo, exhi, &
      	ey, eylo, eyhi, &
      	ez, ezlo, ezhi, &
      	f, flo, fhi, &
      	dtsdx_c2, dtsdy_c2, dtsdz_c2)
    ELSE IF (maxwell_fdtd_solver_id .eq. 1) THEN
      ! Cole-Karkkainen FDTD solver
      CALL WRPX_PXR_PUSH_EVEC_F_CKC( &
        xlo, xhi, ylo, yhi, zlo, zhi, &
      	ex, exlo, exhi, &
      	ey, eylo, eyhi, &
      	ez, ezlo, ezhi, &
      	f, flo, fhi, &
      	dtsdx_c2, dtsdy_c2, dtsdz_c2)
    ENDIF
  end subroutine warpx_push_evec_f

end module warpx_to_pxr_module
