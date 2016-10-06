
module warpx_to_pxr_module

  use iso_c_binding
  use bl_fort_module, only : c_real

  implicit none
  
contains

  subroutine warpx_geteb3d_energy_conserving(np,xp,yp,zp, &
       ex,ey,ez,bx,by,bz,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
       nox,noy,noz,exg,eyg,ezg,bxg,byg,bzg, &
       ll4symtry,l_lower_order_in_v, &  !!!!!!!
       field_gathe_algo) &
       bind(C, name="warpx_geteb3d_energy_conserving")

    interface
       subroutine geteb3d_energy_conserving(np,xp,yp,zp, &
            ex,ey,ez,bx,by,bz,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
            nox,noy,noz,exg,eyg,ezg,bxg,byg,bzg, &
            ll4symtry,l_lower_order_in_v, &
            field_gathe_algo) &
            bind(C, name="geteb3d_energy_conserving")
         use iso_c_binding
         use bl_fort_module, only : c_real
         implicit none
         integer(c_long), intent(in) :: field_gathe_algo
         integer(c_long), intent(in) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
         logical, intent(in)      :: ll4symtry,l_lower_order_in_v
         real(c_real), dimension(np) :: xp,yp,zp,ex,ey,ez,bx,by,bz
         real(c_real), intent(in)    :: xmin,ymin,zmin,dx,dy,dz
         real(c_real), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
         real(c_real), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg
       end subroutine geteb3d_energy_conserving
    end interface

    integer(c_long), intent(in) :: field_gathe_algo
    integer(c_long), intent(in) :: np,nx,ny,nz,nox,noy,noz,nxguard,nyguard,nzguard
    integer(c_int), intent(in)      :: ll4symtry,l_lower_order_in_v
    real(c_real), dimension(np) :: xp,yp,zp,ex,ey,ez,bx,by,bz
    real(c_real), intent(in)    :: xmin,ymin,zmin,dx,dy,dz
    real(c_real), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: exg,eyg,ezg
    real(c_real), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: bxg,byg,bzg

    logical :: pxr_ll4symtry, pxr_l_lower_order_in_v
    
    pxr_ll4symtry = ll4symtry .eq. 1
    pxr_l_lower_order_in_v = l_lower_order_in_v .eq. 1
    
    call geteb3d_energy_conserving(np,xp,yp,zp, &
         ex,ey,ez,bx,by,bz,xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard, &
         nox,noy,noz,exg,eyg,ezg,bxg,byg,bzg, &
         pxr_ll4symtry, pxr_l_lower_order_in_v, &
         field_gathe_algo)

  end subroutine warpx_geteb3d_energy_conserving


  subroutine warpx_pxrpush_em3d_evec_norder(ex,ey,ez,bx,by,bz,jx,jy,jz,mudt,    &
       dtsdx,dtsdy,dtsdz,nx,ny,nz,   &
       norderx,nordery,norderz,             &
       nxguard,nyguard,nzguard,nxs,nys,nzs, &
       l_nodalgrid) &  !!!!!
       bind(C, name="warpx_pxrpush_em3d_evec_norder")
    
    interface
       subroutine pxrpush_em3d_evec_norder(ex,ey,ez,bx,by,bz,jx,jy,jz,mudt, &
            dtsdx,dtsdy,dtsdz,nx,ny,nz,   &
            norderx,nordery,norderz,             &
            nxguard,nyguard,nzguard,nxs,nys,nzs, &
            l_nodalgrid) &  !!!!!
            bind(C, name="pxrpush_em3d_evec_norder")
         use iso_c_binding
         use bl_fort_module, only : c_real
         implicit none
         integer(c_long) :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs,norderx,nordery,norderz
         real(c_real), intent(IN OUT), dimension(-nxguard:nx+nxguard,&
              &                                  -nyguard:ny+nyguard,&
              &                                  -nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
         real(c_real), intent(IN), dimension(-nxguard:nx+nxguard,&
              &                              -nyguard:ny+nyguard,&
              &                              -nzguard:nz+nzguard) :: jx, jy, jz
         real(c_real), intent(IN) :: mudt,dtsdx(norderx/2),dtsdy(nordery/2),dtsdz(norderz/2)
         integer(c_long) :: j,k,l,ist
         logical :: l_nodalgrid
       end subroutine pxrpush_em3d_evec_norder
    end interface

    integer(c_long) :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs,norderx,nordery,norderz
    real(c_real), intent(IN OUT), dimension(-nxguard:nx+nxguard,&
         &                                  -nyguard:ny+nyguard,&
         &                                  -nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
    real(c_real), intent(IN), dimension(-nxguard:nx+nxguard,&
         &                              -nyguard:ny+nyguard,&
         &                              -nzguard:nz+nzguard) :: jx, jy, jz
    real(c_real), intent(IN) :: mudt,dtsdx(norderx/2),dtsdy(nordery/2),dtsdz(norderz/2)
    integer(c_long) :: j,k,l,ist
    integer(c_int) :: l_nodalgrid

    logical pxr_l_nodalgrid

    pxr_l_nodalgrid = l_nodalgrid .eq. 1
    
    call pxrpush_em3d_evec_norder(ex,ey,ez,bx,by,bz,jx,jy,jz,mudt, &
         dtsdx,dtsdy,dtsdz,nx,ny,nz,   &
         norderx,nordery,norderz,             &
         nxguard,nyguard,nzguard,nxs,nys,nzs, &
         pxr_l_nodalgrid)
  end subroutine warpx_pxrpush_em3d_evec_norder


  subroutine warpx_pxrpush_em3d_bvec_norder(ex,ey,ez,bx,by,bz,                  &
       dtsdx,dtsdy,dtsdz,nx,ny,nz,          &
       norderx,nordery,norderz,             &
       nxguard,nyguard,nzguard,nxs,nys,nzs, &
       l_nodalgrid) &  !!!!!
       bind(C, name="warpx_pxrpush_em3d_bvec_norder")

    interface
       subroutine pxrpush_em3d_bvec_norder(ex,ey,ez,bx,by,bz,                  &
            dtsdx,dtsdy,dtsdz,nx,ny,nz,          &
            norderx,nordery,norderz,             &
            nxguard,nyguard,nzguard,nxs,nys,nzs, &
            l_nodalgrid) &  !!!!!
            bind(C, name="pxrpush_em3d_bvec_norder")
         use iso_c_binding
         use bl_fort_module, only : c_real
         implicit none
         integer(c_long) :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs,norderx,nordery,norderz
         real(c_real), intent(IN OUT), dimension(-nxguard:nx+nxguard,&
              &                                  -nyguard:ny+nyguard,&
              &                                  -nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
         real(c_real), intent(IN) :: dtsdx(norderx/2),dtsdy(nordery/2),dtsdz(norderz/2)
         integer(c_long) :: j,k,l,ist
         logical :: l_nodalgrid
       end subroutine pxrpush_em3d_bvec_norder
    end interface

    integer(c_long) :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs,norderx,nordery,norderz
    real(c_real), intent(IN OUT), dimension(-nxguard:nx+nxguard,&
         &                                  -nyguard:ny+nyguard,&
         &                                  -nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
    real(c_real), intent(IN) :: dtsdx(norderx/2),dtsdy(nordery/2),dtsdz(norderz/2)
    integer(c_long) :: j,k,l,ist
    integer(c_int) :: l_nodalgrid

    logical pxr_l_nodalgrid

    pxr_l_nodalgrid = l_nodalgrid .eq. 1

    call pxrpush_em3d_bvec_norder(ex,ey,ez,bx,by,bz,                  &
         dtsdx,dtsdy,dtsdz,nx,ny,nz,          &
         norderx,nordery,norderz,             &
         nxguard,nyguard,nzguard,nxs,nys,nzs, &
         pxr_l_nodalgrid)
  end subroutine warpx_pxrpush_em3d_bvec_norder

end module warpx_to_pxr_module
