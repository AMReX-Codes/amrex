#if (AMREX_SPACEDIM == 3)

#define WRPX_PXR_CURRENT_DEPOSITION       depose_jxjyjz_generic

#elif (AMREX_SPACEDIM == 2)

#ifdef WARPX_DIM_RZ

#define WRPX_PXR_CURRENT_DEPOSITION       depose_jrjtjz_generic_rz

#else

#define WRPX_PXR_CURRENT_DEPOSITION       depose_jxjyjz_generic_2d

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
    l_nodal,lvect,current_depo_algo) &
    bind(C, name="warpx_current_deposition")

    integer, intent(in) :: jx_ntot(AMREX_SPACEDIM), jy_ntot(AMREX_SPACEDIM), jz_ntot(AMREX_SPACEDIM)
    integer(c_long), intent(in) :: jx_ng, jy_ng, jz_ng
    integer(c_long), intent(IN)                                  :: np
    integer(c_long), intent(IN)                                  :: nox,noy,noz
    integer(c_int), intent(in)                                   :: l_nodal
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
    logical(pxr_logical)                                          :: pxr_l_nodal

    ! Compute the number of valid cells and guard cells
    integer(c_long) :: jx_nvalid(AMREX_SPACEDIM), jy_nvalid(AMREX_SPACEDIM), jz_nvalid(AMREX_SPACEDIM), &
                       jx_nguards(AMREX_SPACEDIM), jy_nguards(AMREX_SPACEDIM), jz_nguards(AMREX_SPACEDIM)
    jx_nvalid = jx_ntot - 2*jx_ng
    jy_nvalid = jy_ntot - 2*jy_ng
    jz_nvalid = jz_ntot - 2*jz_ng
    jx_nguards = jx_ng
    jy_nguards = jy_ng
    jz_nguards = jz_ng
    pxr_l_nodal = l_nodal .eq. 1

! Dimension 3
#if (AMREX_SPACEDIM==3)
   CALL WRPX_PXR_CURRENT_DEPOSITION(        &
        jx,jx_nguards,jx_nvalid,            &
        jy,jy_nguards,jy_nvalid,            &
        jz,jz_nguards,jz_nvalid,            &
        np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q, &
        xmin,ymin,zmin,dt,dx,dy,dz,         &
        nox,noy,noz,pxr_l_nodal,current_depo_algo)
! Dimension 2
#elif (AMREX_SPACEDIM==2)
        CALL WRPX_PXR_CURRENT_DEPOSITION(   &
        jx,jx_nguards,jx_nvalid,            &
        jy,jy_nguards,jy_nvalid,            &
        jz,jz_nguards,jz_nvalid,            &
        np,xp,yp,zp,uxp,uyp,uzp,gaminv,w,q, &
        xmin,zmin,dt,dx,dz,nox,noz,pxr_l_nodal, &
        lvect,current_depo_algo)
#endif

  end subroutine warpx_current_deposition

end module warpx_to_pxr_module
