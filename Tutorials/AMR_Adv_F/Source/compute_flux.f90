
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

  subroutine compute_flux(mla,phi,flux,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi(:)
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:), dt(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: dm, ng_p, ng_f, i, n, nlevs

    dm    = mla%dim
    nlevs = mla%nlevel

    ng_p = phi(1)%ng
    ng_f = flux(1,1)%ng

    do n=1,nlevs
       call compute_flux_single_level(mla,phi(n),flux(:,n),dx(n),dt(n),the_bc_tower%bc_tower_array(n))
    end do

    ! set level n-1 fluxes to be the average of the level n fluxes covering it
    ! the loop over nlevs must count backwards to make sure the finer grids are done first
    do n=nlevs,2,-1
       do i=1,dm
          call ml_edge_restriction_c(flux(i,n-1),1,flux(i,n),1,mla%mba%rr(n-1,:),i,1)
       end do
    end do

  end subroutine compute_flux

  subroutine compute_flux_single_level(mla,phi,flux,dx,dt,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi
    type(multifab) , intent(inout) :: flux(:)
    real(kind=dp_t), intent(in   ) :: dx, dt
    type(bc_level) , intent(in   ) :: the_bc_level

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, ng_f, i

    real(kind=dp_t), pointer ::  pp(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)

    dm    = mla%dim

    ng_p = phi%ng
    ng_f = flux(1)%ng

    do i=1,nfabs(phi)
       pp  => dataptr(phi,i)
       fxp => dataptr(flux(1),i)
       fyp => dataptr(flux(2),i)
       lo = lwb(get_box(phi,i))
       hi = upb(get_box(phi,i))
       select case(dm)
       case (2)
          call compute_flux_2d(pp(:,:,1,1), ng_p, &
                               fxp(:,:,1,1),  fyp(:,:,1,1), ng_f, &
                               lo, hi, dx, dt, &
                               the_bc_level%adv_bc_level_array(i,:,:,1))

       case (3)
          fzp => dataptr(flux(3),i)
          call compute_flux_3d(pp(:,:,:,1), ng_p, &
                               fxp(:,:,:,1),  fyp(:,:,:,1), fzp(:,:,:,1), ng_f, &
                               lo, hi, dx, dt, &
                               the_bc_level%adv_bc_level_array(i,:,:,1))
       end select
    end do

  end subroutine compute_flux_single_level

  subroutine compute_flux_2d(phi, ng_p, fluxx, fluxy, ng_f, lo, hi, dx, dt, adv_bc)
    
    use prob_module , only : mu, uadv, vadv
    use slope_module, only : slope_2d

    integer          :: lo(2), hi(2), ng_p, ng_f
    double precision ::   phi(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: dx, dt
    integer          :: adv_bc(:,:)

    ! local variables
    integer          :: i,j
    double precision :: dlft,drgt,dcen,dlim,phi_edge,eps
    double precision, allocatable :: slope(:,:)

    ! ****************************************************************
    ! Advective Flxues
    ! ****************************************************************

    allocate(slope(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))

    eps = 1.d-8 * max(uadv,vadv)

    ! Compute the slopes in the x-direction
    call slope_2d(phi  ,lo(1)-ng_p,lo(2)-ng_p,hi(1)+ng_p,hi(2)+ng_p, &
                  slope,lo(1)-1   ,lo(2)-1   ,hi(1)+1   ,hi(2)+1, &
                  lo(1),lo(2),hi(1),hi(2),1,1) 

    ! x-fluxes: add the advective fluxes -- these are 2nd order slopes
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1

          if (abs(uadv) < eps) then
              phi_edge = 0.5d0 * (phi(i,j) + phi(i-1,j))
          else if (uadv > 0.d0) then
              phi_edge = phi(i-1,j) + 0.5d0 * (1.d0 - uadv*dt/dx)*slope(i-1,j)
          else
              phi_edge = phi(i  ,j) - 0.5d0 * (1.d0 + uadv*dt/dx)*slope(i  ,j)
          endif

          fluxx(i,j) = fluxx(i,j) - uadv * phi_edge

       end do
    end do

    ! Compute the slopes in the y-direction
    call slope_2d(phi  ,lo(1)-ng_p,lo(2)-ng_p,hi(1)+ng_p,hi(2)+ng_p, &
                  slope,lo(1)-1   ,lo(2)-1   ,hi(1)+1   ,hi(2)+1, &
                  lo(1),lo(2),hi(1),hi(2),1,2) 

    ! y-fluxes: add the advective fluxes
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          
          if (abs(vadv) < eps) then
              phi_edge = 0.5d0 * (phi(i,j) + phi(i,j-1))
          else if (vadv > 0.d0) then
              phi_edge = phi(i,j-1) + 0.5d0 * (1.d0 - vadv*dt/dx)*slope(i,j-1)
          else
              phi_edge = phi(i,j  ) - 0.5d0 * (1.d0 + vadv*dt/dx)*slope(i,j  )
          endif

          fluxy(i,j) = fluxy(i,j) - vadv * phi_edge

       end do
    end do

    deallocate(slope)

  end subroutine compute_flux_2d

  subroutine compute_flux_3d(phi, ng_p, fluxx, fluxy, fluxz, ng_f, &
                             lo, hi, dx, dt, adv_bc)

    use prob_module , only : mu, uadv, vadv, wadv
    use slope_module, only : slope_3d

    integer          :: lo(3), hi(3), ng_p, ng_f
    double precision ::   phi(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: dx, dt
    integer          :: adv_bc(:,:)

    ! local variables
    integer          :: i,j,k
    double precision :: dlft,drgt,dcen,dlim,phi_edge,eps
    double precision, allocatable :: slope(:,:,:)

    ! ****************************************************************
    ! Advective fluxes
    ! ****************************************************************

    allocate(slope(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    eps = 1.d-8 * max(max(uadv,vadv),wadv)

    ! Compute the slopes in the x-direction
    call slope_3d(phi  ,lo(1)-ng_p,lo(2)-ng_p,lo(3)-ng_p,hi(1)+ng_p,hi(2)+ng_p,hi(3)+ng_p,&
                  slope,lo(1)-1   ,lo(2)-1   ,lo(3)-1   ,hi(1)+1   ,hi(2)+1   ,hi(3)+1,   &
                  lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),1,1) 

    ! x-fluxes: add the advective fluxes -- these are 2nd order slopes
    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1

             if (abs(uadv) < eps) then
                 phi_edge = 0.5d0 * (phi(i,j,k) + phi(i-1,j,k))
             else if (uadv > 0.d0) then
                 phi_edge = phi(i-1,j,k) + 0.5d0 * (1.d0 - uadv*dt/dx)*slope(i-1,j,k)
             else
                 phi_edge = phi(i  ,j,k) - 0.5d0 * (1.d0 + uadv*dt/dx)*slope(i  ,j,k)
             endif
   
             fluxx(i,j,k) = fluxx(i,j,k) - uadv * phi_edge

          end do
       end do
    end do
    !$omp end parallel do

    ! Compute the slopes in the y-direction
    call slope_3d(phi  ,lo(1)-ng_p,lo(2)-ng_p,lo(3)-ng_p,hi(1)+ng_p,hi(2)+ng_p,hi(3)+ng_p,&
                  slope,lo(1)-1   ,lo(2)-1   ,lo(3)-1   ,hi(1)+1   ,hi(2)+1   ,hi(3)+1,   &
                  lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),1,2) 

    ! y-fluxes: add the advective fluxes
    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)

             if (abs(vadv) < eps) then
                 phi_edge = 0.5d0 * (phi(i,j,k) + phi(i,j-1,k))
             else if (vadv > 0.d0) then
                 phi_edge = phi(i,j-1,k) + 0.5d0 * (1.d0 - vadv*dt/dx)*slope(i,j-1,k)
             else
                 phi_edge = phi(i,j  ,k) - 0.5d0 * (1.d0 + vadv*dt/dx)*slope(i,j  ,k)
             endif
   
             fluxy(i,j,k) = fluxy(i,j,k) - vadv * phi_edge
   
          end do
       end do
    end do
    !$omp end parallel do

    ! Compute the slopes in the y-direction
    call slope_3d(phi  ,lo(1)-ng_p,lo(2)-ng_p,lo(3)-ng_p,hi(1)+ng_p,hi(2)+ng_p,hi(3)+ng_p,&
                  slope,lo(1)-1   ,lo(2)-1   ,lo(3)-1   ,hi(1)+1   ,hi(2)+1   ,hi(3)+1,   &
                  lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),1,3) 

    ! z-fluxes: add the advective fluxes
    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             if (abs(wadv) < eps) then
                 phi_edge = 0.5d0 * (phi(i,j,k) + phi(i,j,k-1))
             else if (wadv > 0.d0) then
                 phi_edge = phi(i,j,k-1) + 0.5d0 * (1.d0 - wadv*dt/dx)*slope(i,j,k-1)
             else
                 phi_edge = phi(i,j,k  ) - 0.5d0 * (1.d0 + wadv*dt/dx)*slope(i,j,k  )
             endif
   
             fluxz(i,j,k) = fluxz(i,j,k) - wadv * phi_edge
   
          end do
       end do
    end do
    !$omp end parallel do

    deallocate(slope)

  end subroutine compute_flux_3d

end module compute_flux_module

