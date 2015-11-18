
module compute_flux_module

  use multifab_module
  use ml_layout_module
  use multifab_physbc_module
  use multifab_fill_ghost_module
  use ml_cc_restriction_module
  use bndry_reg_module
  use define_bc_module
  use bc_module

  implicit none

  private

  public :: compute_flux, compute_flux_single_level

contains

  subroutine compute_flux(mla,phi,velocity,flux,dx,dt)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi(:)
    type(multifab) , intent(in   ) :: velocity(:,:)
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:), dt(:)

    ! local variables
    integer :: dm, ng_p, ng_f, i, n, nlevs

    dm    = mla%dim
    nlevs = mla%nlevel

    ng_p = phi(1)%ng
    ng_f = flux(1,1)%ng

    do n=1,nlevs
       call compute_flux_single_level(mla,phi(n),velocity(n,:),flux(n,:),dx(n),dt(n))
    end do

    ! set level n-1 fluxes to be the average of the level n fluxes covering it
    ! the loop over nlevs must count backwards to make sure the finer grids are done first
    do n=nlevs,2,-1
       do i=1,dm
          call ml_edge_restriction_c(flux(n-1,i),1,flux(n,i),1,mla%mba%rr(n-1,:),i,1)
       end do
    end do

  end subroutine compute_flux

  subroutine compute_flux_single_level(mla,phi,velocity,flux,dx,dt)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi
    type(multifab) , intent(in   ) :: velocity(:)
    type(multifab) , intent(inout) :: flux(:)
    real(kind=dp_t), intent(in   ) :: dx, dt

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, ng_f, ng_u, i

    real(kind=dp_t), pointer ::  pp(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)
    real(kind=dp_t), pointer ::  up(:,:,:,:)
    real(kind=dp_t), pointer ::  vp(:,:,:,:)
    real(kind=dp_t), pointer ::  wp(:,:,:,:)

    dm    = mla%dim

    ng_p = phi%ng
    ng_f = flux(1)%ng
    ng_u = velocity(1)%ng

    do i=1,nfabs(phi)
       pp  => dataptr(phi,i)
       fxp => dataptr(flux(1),i)
       fyp => dataptr(flux(2),i)
       up  => dataptr(velocity(1),i)
       vp  => dataptr(velocity(2),i)
       lo = lwb(get_box(phi,i))
       hi = upb(get_box(phi,i))
       select case(dm)
       case (2)
          call compute_flux_2d(pp(:,:,1,1), ng_p, &
                               up(:,:,1,1), vp(:,:,1,1), ng_u, &
                               fxp(:,:,1,1),  fyp(:,:,1,1), ng_f, &
                               lo, hi, dx, dt)

       case (3)
          fzp => dataptr(flux(3),i)
          wp  => dataptr(velocity(3),i)
          call compute_flux_3d(pp(:,:,:,1), ng_p, &
                               up(:,:,:,1), vp(:,:,:,1), wp(:,:,:,1), ng_u, &
                               fxp(:,:,:,1),  fyp(:,:,:,1), fzp(:,:,:,1), ng_f, &
                               lo, hi, dx, dt)
       end select
    end do

  end subroutine compute_flux_single_level

  subroutine compute_flux_2d(phi, ng_p, umac, vmac, ng_u, fluxx, fluxy, ng_f, &
                             lo, hi, dx, dt)
    
    use slope_module, only : slope_2d

    integer          :: lo(2), hi(2), ng_p, ng_f, ng_u
    double precision ::   phi(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision ::  umac(lo(1)-ng_u:,lo(2)-ng_u:)
    double precision ::  vmac(lo(1)-ng_u:,lo(2)-ng_u:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: dx, dt

    ! local variables
    integer          :: i,j
    double precision :: hdtdx

    double precision, allocatable :: slope(:,:)
    double precision, allocatable :: phix   (:,:)
    double precision, allocatable :: phiy   (:,:)

    allocate(slope  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))

    ! normal (1D) predictor states
    ! allocated from lo  :hi+1 in the normal direction
    !                lo-1:hi+1 in the transverse directions
    allocate(phix(lo(1)  :hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(phiy(lo(1)-1:hi(1)+1,lo(2)  :hi(2)+1))

    hdtdx = dt/(2.d0*dx)

    ! ****************************************************************
    ! Advective Flxues
    ! ****************************************************************

    ! Compute the slopes in the x-direction
    call slope_2d(phi  ,lo(1)-ng_p,lo(2)-ng_p,hi(1)+ng_p,hi(2)+ng_p, &
                  slope,lo(1)-1   ,lo(2)-1   ,hi(1)+1   ,hi(2)+1, &
                  lo(1),lo(2),hi(1),hi(2),1,1) 

    ! compute phi on x faces using umac to upwind; ignore transverse terms
    do j=lo(2)-1,hi(2)+1
       do i=lo(1),hi(1)+1

          if (umac(i,j) .lt. 0.d0) then
             phix(i,j) = phi(i  ,j) - (0.5d0 + hdtdx*umac(i,j))*slope(i  ,j)
          else
             phix(i,j) = phi(i-1,j) + (0.5d0 - hdtdx*umac(i,j))*slope(i-1,j)
          end if

       end do
    end do

    ! compute the slopes in the y-direction
    call slope_2d(phi  ,lo(1)-ng_p,lo(2)-ng_p,hi(1)+ng_p,hi(2)+ng_p, &
                  slope,lo(1)-1   ,lo(2)-1   ,hi(1)+1   ,hi(2)+1, &
                  lo(1),lo(2),hi(1),hi(2),1,2) 

    ! compute phi on y faces using umac to upwind; ignore transverse terms
    do j=lo(2),hi(2)+1
       do i=lo(1)-1,hi(1)+1

          if (vmac(i,j) .lt. 0.d0) then
             phiy(i,j) = phi(i,j  ) - (0.5d0 + hdtdx*vmac(i,j))*slope(i,j  )
          else
             phiy(i,j) = phi(i,j-1) + (0.5d0 - hdtdx*vmac(i,j))*slope(i,j-1)
          end if

       end do
    end do

    ! update phi on x faces by adding in y-transverse terms
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1

          if (umac(i,j) .lt. 0.d0) then
             fluxx(i,j) = phix(i,j) &
                  - hdtdx*( 0.5d0*(vmac(i  ,j+1)+vmac(i  ,j)) * (phiy(i  ,j+1)-phiy(i  ,j)) )
          else
             fluxx(i,j) = phix(i,j) &
                  - hdtdx*( 0.5d0*(vmac(i-1,j+1)+vmac(i-1,j)) * (phiy(i-1,j+1)-phiy(i-1,j)) )
          end if

          ! compute final x-fluxes
          fluxx(i,j) = fluxx(i,j)*umac(i,j)

       end do
    end do

    ! update phi on y faces by adding in x-transverse terms
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)

          if (vmac(i,j) .lt. 0.d0) then
             fluxy(i,j) = phiy(i,j) &
                  - hdtdx*( 0.5d0*(umac(i+1,j  )+umac(i,j  )) * (phix(i+1,j  )-phix(i,j  )) )
          else
             fluxy(i,j) = phiy(i,j) &
                  - hdtdx*( 0.5d0*(umac(i+1,j-1)+umac(i,j-1)) * (phix(i+1,j-1)-phix(i,j-1)) )
          end if

          ! compute final y-fluxes
          fluxy(i,j) = fluxy(i,j)*vmac(i,j)

       end do
    end do

    deallocate(slope,phix,phiy)

  end subroutine compute_flux_2d

  subroutine compute_flux_3d(phi, ng_p, umac, vmac, wmac, ng_u, fluxx, fluxy, fluxz, ng_f, &
                             lo, hi, dx, dt)

    use slope_module, only : slope_3d

    integer          :: lo(3), hi(3), ng_p, ng_f, ng_u
    double precision ::   phi(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision ::  umac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    double precision ::  vmac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    double precision ::  wmac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: dx, dt

    ! local variables
    integer          :: i,j,k
    double precision :: hdtdx, tdtdx
    double precision, allocatable :: slope(:,:,:)
    double precision, allocatable :: phix   (:,:,:)
    double precision, allocatable :: phix_y (:,:,:)
    double precision, allocatable :: phix_z (:,:,:)
    double precision, allocatable :: phiy   (:,:,:)
    double precision, allocatable :: phiy_x (:,:,:)
    double precision, allocatable :: phiy_z (:,:,:)
    double precision, allocatable :: phiz   (:,:,:)
    double precision, allocatable :: phiz_x (:,:,:)
    double precision, allocatable :: phiz_y (:,:,:)

    allocate(slope  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    ! normal (1D) predictor states
    ! allocated from lo  :hi+1 in the normal direction
    !                lo-1:hi+1 in the transverse directions
    allocate(phix(lo(1)  :hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(phiy(lo(1)-1:hi(1)+1,lo(2)  :hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(phiz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)  :hi(3)+1))

    ! These are transverse terms.  The size allocation is tricky.
    ! lo  :hi+1 in the normal direction
    ! lo  :hi   in the transverse direction
    ! lo-1:hi+1 in the unused direction
    allocate(phix_y (lo(1)  :hi(1)+1,lo(2)  :hi(2)  ,lo(3)-1:hi(3)+1))
    allocate(phix_z (lo(1)  :hi(1)+1,lo(2)-1:hi(2)+1,lo(3)  :hi(3)  ))
    allocate(phiy_x (lo(1)  :hi(1)  ,lo(2)  :hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(phiy_z (lo(1)-1:hi(1)+1,lo(2)  :hi(2)+1,lo(3)  :hi(3)  ))
    allocate(phiz_x (lo(1)  :hi(1)  ,lo(2)-1:hi(2)+1,lo(3)  :hi(3)+1))
    allocate(phiz_y (lo(1)-1:hi(1)+1,lo(2)  :hi(2)  ,lo(3)  :hi(3)+1))

    hdtdx = dt/(2.d0*dx)
    tdtdx = dt/(3.d0*dx)

    ! Compute the slopes in the x-direction
    call slope_3d(phi  ,lo(1)-ng_p,lo(2)-ng_p,lo(3)-ng_p,hi(1)+ng_p,hi(2)+ng_p,hi(3)+ng_p, &
                  slope,lo(1)-1   ,lo(2)-1   ,lo(3)-1   ,hi(1)+1   ,hi(2)+1   ,hi(3)+1, &
                  lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),1,1) 
    
    ! compute phi on x faces using umac to upwind; ignore transverse terms
    do k=lo(3)-1,hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1),hi(1)+1

             if (umac(i,j,k) .lt. 0.d0) then
                phix(i,j,k) = phi(i  ,j,k) - (0.5d0 + hdtdx*umac(i,j,k))*slope(i  ,j,k)
             else
                phix(i,j,k) = phi(i-1,j,k) + (0.5d0 - hdtdx*umac(i,j,k))*slope(i-1,j,k)
             end if

          end do
       end do
    end do

    call slope_3d(phi  ,lo(1)-ng_p,lo(2)-ng_p,lo(3)-ng_p,hi(1)+ng_p,hi(2)+ng_p,hi(3)+ng_p, &
                  slope,lo(1)-1   ,lo(2)-1   ,lo(3)-1   ,hi(1)+1   ,hi(2)+1   ,hi(3)+1, &
                  lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),1,2) 

    ! compute phi on y faces using umac to upwind; ignore transverse terms
    do k=lo(3)-1,hi(3)+1
       do j=lo(2),hi(2)+1
          do i=lo(1)-1,hi(1)+1

             if (vmac(i,j,k) .lt. 0.d0) then
                phiy(i,j,k) = phi(i,j  ,k) - (0.5d0 + hdtdx*vmac(i,j,k))*slope(i,j  ,k)
             else
                phiy(i,j,k) = phi(i,j-1,k) + (0.5d0 - hdtdx*vmac(i,j,k))*slope(i,j-1,k)
             end if

          end do
       end do
    end do

    call slope_3d(phi  ,lo(1)-ng_p,lo(2)-ng_p,lo(3)-ng_p,hi(1)+ng_p,hi(2)+ng_p,hi(3)+ng_p, &
                  slope,lo(1)-1   ,lo(2)-1   ,lo(3)-1   ,hi(1)+1   ,hi(2)+1   ,hi(3)+1, &
                  lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),1,3) 

    ! compute phi on z faces using umac to upwind; ignore transverse terms
    do k=lo(3),hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1

             if (wmac(i,j,k) .lt. 0.d0) then
                phiz(i,j,k) = phi(i,j,k  ) - (0.5d0 + hdtdx*wmac(i,j,k))*slope(i,j,k  )
             else
                phiz(i,j,k) = phi(i,j,k-1) + (0.5d0 - hdtdx*wmac(i,j,k))*slope(i,j,k-1)
             end if

          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!
    ! transverse terms
    !!!!!!!!!!!!!!!!!!!!

    ! update phi on x faces by adding in y-transverse terms
    do k=lo(3)-1,hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1

             if (umac(i,j,k) .lt. 0.d0) then
                phix_y(i,j,k) = phix(i,j,k) &
                     - tdtdx* (0.5d0*(vmac(i  ,j+1,k)+vmac(i  ,j,k)) * (phiy(i  ,j+1,k)-phiy(i  ,j,k)) )
             else
                phix_y(i,j,k) = phix(i,j,k) &
                     - tdtdx* (0.5d0*(vmac(i-1,j+1,k)+vmac(i-1,j,k)) * (phiy(i-1,j+1,k)-phiy(i-1,j,k)) )
             end if

          end do
       end do
    end do

    ! update phi on x faces by adding in z-transverse terms
    do k=lo(3),hi(3)
       do j=lo(2)-1,hi(2)+1
          do i=lo(1),hi(1)+1

             if (umac(i,j,k) .lt. 0.d0) then
                phix_z(i,j,k) = phix(i,j,k) &
                     - tdtdx* (0.5d0*(wmac(i  ,j,k+1)+wmac(i  ,j,k)) * (phiz(i  ,j,k+1)-phiz(i  ,j,k)) )
             else
                phix_z(i,j,k) = phix(i,j,k) &
                     - tdtdx* (0.5d0*(wmac(i-1,j,k+1)+wmac(i-1,j,k)) * (phiz(i-1,j,k+1)-phiz(i-1,j,k)) )
             end if

          end do
       end do
    end do

    ! update phi on y faces by adding in x-transverse terms
    do k=lo(3)-1,hi(3)+1
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)

             if (vmac(i,j,k) .lt. 0.d0) then
                phiy_x(i,j,k) = phiy(i,j,k) &
                     - tdtdx* (0.5d0*(umac(i+1,j  ,k)+umac(i,j  ,k)) * (phix(i+1,j  ,k)-phix(i,j  ,k)) )
             else
                phiy_x(i,j,k) = phiy(i,j,k) &
                     - tdtdx* (0.5d0*(umac(i+1,j-1,k)+umac(i,j-1,k)) * (phix(i+1,j-1,k)-phix(i,j-1,k)) )
             end if

          end do
       end do
    end do

    ! update phi on y faces by adding in z-transverse terms
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1)-1,hi(1)+1

             if (vmac(i,j,k) .lt. 0.d0) then
                phiy_z(i,j,k) = phiy(i,j,k) &
                     - tdtdx* (0.5d0*(wmac(i,j  ,k+1)+wmac(i,j  ,k)) * (phiz(i,j  ,k+1)-phiz(i,j  ,k)) )
             else
                phiy_z(i,j,k) = phiy(i,j,k) &
                     - tdtdx* (0.5d0*(wmac(i,j-1,k+1)+wmac(i,j-1,k)) * (phiz(i,j-1,k+1)-phiz(i,j-1,k)) )
             end if

          end do
       end do
    end do

    ! update phi on z faces by adding in x-transverse terms
    do k=lo(3),hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1),hi(1)

             if (wmac(i,j,k) .lt. 0.d0) then
                phiz_x(i,j,k) = phiz(i,j,k) &
                     - tdtdx* (0.5d0*(umac(i+1,j,k  )+umac(i,j,k  )) * (phix(i+1,j,k  )-phix(i,j,k  )) )
             else
                phiz_x(i,j,k) = phiz(i,j,k) &
                     - tdtdx* (0.5d0*(umac(i+1,j,k-1)+umac(i,j,k-1)) * (phix(i+1,j,k-1)-phix(i,j,k-1)) )
             end if

          end do
       end do
    end do

    ! update phi on z faces by adding in y-transverse terms
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1)-1,hi(1)+1

             if (wmac(i,j,k) .lt. 0.d0) then
                phiz_y(i,j,k) = phiz(i,j,k) &
                     - tdtdx* (0.5d0*(vmac(i,j+1,k  )+vmac(i,j,k  )) * (phiy(i,j+1,k  )-phiy(i,j,k  )) )
             else
                phiz_y(i,j,k) = phiz(i,j,k) &
                     - tdtdx* (0.5d0*(vmac(i,j+1,k-1)+vmac(i,j,k-1)) * (phiy(i,j+1,k-1)-phiy(i,j,k-1)) )
             end if

          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!
    ! final edge states
    !!!!!!!!!!!!!!!!!!!!

    ! update phi on x faces by adding in yz and zy transverse terms
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1

             if (umac(i,j,k) .lt. 0.d0) then
                fluxx(i,j,k) = phix(i,j,k) &
                     - hdtdx*( 0.5d0*(vmac(i  ,j+1,k  )+vmac(i  ,j,k)) * (phiy_z(i  ,j+1,k  )-phiy_z(i  ,j,k)) ) &
                     - hdtdx*( 0.5d0*(wmac(i  ,j  ,k+1)+wmac(i  ,j,k)) * (phiz_y(i  ,j  ,k+1)-phiz_y(i  ,j,k)) )
             else
                fluxx(i,j,k) = phix(i,j,k) &
                     - hdtdx*( 0.5d0*(vmac(i-1,j+1,k  )+vmac(i-1,j,k)) * (phiy_z(i-1,j+1,k  )-phiy_z(i-1,j,k)) ) &
                     - hdtdx*( 0.5d0*(wmac(i-1,j  ,k+1)+wmac(i-1,j,k)) * (phiz_y(i-1,j  ,k+1)-phiz_y(i-1,j,k)) )
             end if

             ! compute final x-fluxes
             fluxx(i,j,k) = fluxx(i,j,k)*umac(i,j,k)

          end do
       end do
    end do

    ! update phi on y faces by adding in xz and zx transverse terms
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)

             if (vmac(i,j,k) .lt. 0.d0) then
                fluxy(i,j,k) = phiy(i,j,k) &
                     - hdtdx*( 0.5d0*(umac(i+1,j  ,k  )+umac(i,j  ,k)) * (phix_z(i+1,j  ,k  )-phix_z(i,j  ,k)) ) &
                     - hdtdx*( 0.5d0*(wmac(i  ,j  ,k+1)+wmac(i,j  ,k)) * (phiz_x(i  ,j  ,k+1)-phiz_x(i,j  ,k)) )
             else
                fluxy(i,j,k) = phiy(i,j,k) &
                     - hdtdx*( 0.5d0*(umac(i+1,j-1,k  )+umac(i,j-1,k)) * (phix_z(i+1,j-1,k  )-phix_z(i,j-1,k)) ) &
                     - hdtdx*( 0.5d0*(wmac(i  ,j-1,k+1)+wmac(i,j-1,k)) * (phiz_x(i  ,j-1,k+1)-phiz_x(i,j-1,k)) )
             end if

             ! compute final y-fluxes
             fluxy(i,j,k) = fluxy(i,j,k)*vmac(i,j,k)

          end do
       end do
    end do

    ! update phi on z faces by adding in xy and yx transverse terms
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             if (wmac(i,j,k) .lt. 0.d0) then
                fluxz(i,j,k) = phiz(i,j,k) &
                     - hdtdx*( 0.5d0*(umac(i+1,j  ,k  )+umac(i  ,j,k)) * (phix_y(i+1,j  ,k  )-phix_y(i,j,k  )) ) &
                     - hdtdx*( 0.5d0*(vmac(i  ,j+1,k  )+vmac(i  ,j,k)) * (phiy_x(i  ,j+1,k  )-phiy_x(i,j,k  )) )
             else
                fluxz(i,j,k) = phiz(i,j,k) &
                     - hdtdx*( 0.5d0*(umac(i+1,j  ,k-1)+umac(i,j,k-1)) * (phix_y(i+1,j  ,k-1)-phix_y(i,j,k-1)) ) &
                     - hdtdx*( 0.5d0*(vmac(i  ,j+1,k-1)+vmac(i,j,k-1)) * (phiy_x(i  ,j+1,k-1)-phiy_x(i,j,k-1)) )
             end if

             ! compute final z-fluxes
             fluxz(i,j,k) = fluxz(i,j,k)*wmac(i,j,k)

          end do
       end do
    end do


    deallocate(slope,phix,phix_y,phix_z,phiy,phiy_x,phiy_z)
    deallocate(phiz,phiz_x,phiz_y)

  end subroutine compute_flux_3d

end module compute_flux_module

