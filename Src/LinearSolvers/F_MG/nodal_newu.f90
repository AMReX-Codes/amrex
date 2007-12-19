module nodal_newu_module
 
!  use stencil_module
!  use bndry_reg_module
!  use mg_module
!  use ml_boxarray_module
!  use ml_layout_module
!  use bl_mem_stat_module
!  use bl_timer_module
!  use bl_IO_module

  use multifab_module
 
!  use ml_restriction_module
!  use ml_prolongation_module
!  use ml_interface_stencil_module
!  use ml_util_module
 
  implicit none
 
contains

!   ********************************************************************************************* !

    subroutine mkunew(unew,phi,coeff,dx,ng)

      type(multifab) , intent(inout) :: unew
      type(multifab) , intent(in   ) :: phi
      type(multifab) , intent(in   ) :: coeff
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: ng
      integer :: i,dm
 
      real(kind=dp_t), pointer :: up(:,:,:,:) 
      real(kind=dp_t), pointer :: pp(:,:,:,:) 
      real(kind=dp_t), pointer :: rp(:,:,:,:) 

      dm = unew%dim

      do i = 1, unew%nboxes
         if ( multifab_remote(unew, i) ) cycle
         up => dataptr(unew, i)
         pp => dataptr(phi, i)
         rp  => dataptr(coeff, i)
         select case (dm)
            case (2)
              call mkunew_2d(up(:,:,1,:), pp(:,:,1,1) ,rp(:,:,1,1), dx, ng)
            case (3)
              call mkunew_3d(up(:,:,:,:), pp(:,:,:,1), rp(:,:,:,1), dx, ng)
         end select
      end do
      call multifab_fill_boundary(unew)

    end subroutine mkunew

!   ********************************************************************************************* !

    subroutine mkunew_2d(unew,phi,coeff,dx,ng)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) :: unew(-ng:,-ng:,:)
      real(kind=dp_t), intent(inout) :: phi(-1:,-1:)
      real(kind=dp_t), intent(in   ) :: coeff(-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer :: i,j,nx,ny
      real(kind=dp_t) :: gpx,gpy

      nx = size(phi,dim=1)-3
      ny = size(phi,dim=2)-3

      do j = 0,ny-1
      do i = 0,nx-1
         gpx = HALF*(phi(i+1,j) + phi(i+1,j+1) - &
                     phi(i  ,j) - phi(i  ,j+1) ) /dx(1)
         gpy = HALF*(phi(i,j+1) + phi(i+1,j+1) - &
                     phi(i,j  ) - phi(i+1,j  ) ) /dx(2)
         unew(i,j,1) = unew(i,j,1) - gpx * coeff(i,j)
         unew(i,j,2) = unew(i,j,2) - gpy * coeff(i,j)
      end do
      end do

    end subroutine mkunew_2d

!   ********************************************************************************************* !

    subroutine mkunew_3d(unew,phi,rhohalf,dx,ng)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) :: unew(-ng:,-ng:,-ng:,:)
      real(kind=dp_t), intent(inout) :: phi(-ng:,-ng:,-ng:)
      real(kind=dp_t), intent(in   ) :: rhohalf(-1:,-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: nx,ny,nz,i,j,k
      real(kind=dp_t) :: gpx,gpy,gpz

      nx = size(phi,dim=1)-3
      ny = size(phi,dim=2)-3
      nz = size(phi,dim=3)-3

      do k = 0,nz-1
      do j = 0,ny-1
      do i = 0,nx-1
         gpx = FOURTH*(phi(i+1,j,k  ) + phi(i+1,j+1,k  ) &
                      +phi(i+1,j,k+1) + phi(i+1,j+1,k+1) & 
                      -phi(i  ,j,k  ) - phi(i  ,j+1,k  ) &
                      -phi(i  ,j,k+1) - phi(i  ,j+1,k+1) ) /dx(1)
         gpy = FOURTH*(phi(i,j+1,k  ) + phi(i+1,j+1,k  ) &
                      +phi(i,j+1,k+1) + phi(i+1,j+1,k+1) & 
                      -phi(i,j  ,k  ) - phi(i+1,j  ,k  ) &
                      -phi(i,j  ,k+1) - phi(i+1,j  ,k+1) ) /dx(2)
         gpz = FOURTH*(phi(i,j  ,k+1) + phi(i+1,j  ,k+1) &
                      +phi(i,j+1,k+1) + phi(i+1,j+1,k+1) & 
                      -phi(i,j  ,k  ) - phi(i+1,j  ,k  ) &
                      -phi(i,j+1,k  ) - phi(i+1,j+1,k  ) ) /dx(3)

         unew(i,j,k,1) = unew(i,j,k,1) - gpx / rhohalf(i,j,k)
         unew(i,j,k,2) = unew(i,j,k,2) - gpy / rhohalf(i,j,k)
         unew(i,j,k,3) = unew(i,j,k,3) - gpz / rhohalf(i,j,k)
      end do
      end do
      end do

    end subroutine mkunew_3d

end module nodal_newu_module
