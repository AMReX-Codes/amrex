
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
    
    use prob_module, only : mu, uadv, vadv

    integer          :: lo(2), hi(2), ng_p, ng_f
    double precision ::   phi(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: dx, dt
    integer          :: adv_bc(:,:)

    ! local variables
    integer          :: i,j
    double precision :: dlft,drgt,dcen,dlim,slope,phi_edge

    ! ****************************************************************
    ! Initialize with the diffusive fluxes
    ! ****************************************************************

    ! x-fluxes
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          fluxx(i,j) = mu * ( phi(i,j) - phi(i-1,j) ) / dx
       end do
    end do

    ! lo-x diffusive boundary conditions
    if (adv_bc(1,1) .eq. EXT_DIR) then
       i=lo(1)
       do j=lo(2),hi(2)
          ! divide by 0.5*dx since the ghost cell value represents
          ! the value at the wall, not the ghost cell-center
          fluxx(i,j) = mu * ( phi(i,j) - phi(i-1,j) ) / (0.5d0*dx)
       end do
    else if (adv_bc(1,1) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxx(lo(1),lo(2):hi(2)) = 0.d0
    end if

    ! hi-x diffusive boundary conditions
    if (adv_bc(1,2) .eq. EXT_DIR) then
       i=hi(1)+1
       do j=lo(2),hi(2)
          ! divide by 0.5*dx since the ghost cell value represents
          ! the value at the wall, not the ghost cell-center
          fluxx(i,j) = mu * ( phi(i,j) - phi(i-1,j) ) / (0.5d0*dx)
       end do
    else if (adv_bc(1,2) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxx(hi(1)+1,lo(2):hi(2)) = 0.d0
    end if

    ! y-fluxes
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          fluxy(i,j) = mu * ( phi(i,j) - phi(i,j-1) ) / dx
       end do
    end do

    ! lo-y boundary conditions
    if (adv_bc(2,1) .eq. EXT_DIR) then
       j=lo(2)
       do i=lo(1),hi(1)
          ! divide by 0.5*dx since the ghost cell value represents
          ! the value at the wall, not the ghost cell-center
          fluxy(i,j) = mu * ( phi(i,j) - phi(i,j-1) ) / (0.5d0*dx)
       end do
    else if (adv_bc(2,1) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxy(lo(1):hi(1),lo(2)) = 0.d0
    end if

    ! hi-y boundary conditions
    if (adv_bc(2,2) .eq. EXT_DIR) then
       j=hi(2)+1
       do i=lo(1),hi(1)
          ! divide by 0.5*dx since the ghost cell value represents
          ! the value at the wall, not the ghost cell-center
          fluxy(i,j) = mu * ( phi(i,j) - phi(i,j-1) ) / (0.5d0*dx)
       end do
    else if (adv_bc(2,2) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxy(lo(1):hi(1),hi(2)+1) = 0.d0
    end if

    ! ****************************************************************
    ! Add the advective fluxes (here we assume constant velocity so no bc's)
    ! ****************************************************************

    ! x-fluxes: add the advective fluxes -- these are 2nd order slopes
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1

          ! Assume that uadv > , which means we need the slope in cell (i-1,j)
          dlft = 2.d0*(phi(i-1,j) - phi(i-2,j))
          drgt = 2.d0*(phi(i  ,j) - phi(i-1,j))
          dcen = .25d0 * (dlft+drgt)
          if (dlft*drgt .ge. 0.d0) then
             dlim = min( abs(dlft), abs(drgt) )
          else
             dlim = 0.d0
          endif
          slope = sign(1.d0, dcen) * min( dlim, abs(dcen) )
          
          phi_edge = phi(i-1,j) + 0.5d0 * (1.d0 - uadv*dt/dx)*slope

          fluxx(i,j) = fluxx(i,j) - uadv * phi_edge

       end do
    end do

    ! y-fluxes: add the advective fluxes
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)

          ! Assume that vadv > 0, which means we need the slope in (i,j-1)
          dlft = 2.d0*(phi(i,j-1) - phi(i,j-2))
          drgt = 2.d0*(phi(i,j  ) - phi(i,j-1))
          dcen = .25d0 * (dlft+drgt)
          if (dlft*drgt .ge. 0.d0) then
             dlim = min( abs(dlft), abs(drgt) )
          else
             dlim = 0.d0
          endif
          slope = sign(1.d0, dcen) * min( dlim, abs(dcen) )
          
          phi_edge = phi(i,j-1) + 0.5d0 * (1.d0 - vadv*dt/dx) * slope

          fluxy(i,j) = fluxy(i,j) - vadv * phi_edge

       end do
    end do

  end subroutine compute_flux_2d

  subroutine compute_flux_3d(phi, ng_p, fluxx, fluxy, fluxz, ng_f, &
                             lo, hi, dx, dt, adv_bc)

    use prob_module, only : mu, uadv, vadv, wadv

    integer          :: lo(3), hi(3), ng_p, ng_f
    double precision ::   phi(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: dx, dt
    integer          :: adv_bc(:,:)

    ! local variables
    integer          :: i,j,k
    double precision :: dlft,drgt,dcen,dlim,slope,phi_edge

    ! ****************************************************************
    ! Initialize with the diffusive fluxes
    ! ****************************************************************

    ! x-fluxes
    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             fluxx(i,j,k) = mu * ( phi(i,j,k) - phi(i-1,j,k) ) / dx
          end do
       end do
    end do
    !$omp end parallel do

    ! lo-x boundary conditions
    if (adv_bc(1,1) .eq. EXT_DIR) then
       i=lo(1)
       !$omp parallel do private(j,k)
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             ! divide by 0.5*dx since the ghost cell value represents
             ! the value at the wall, not the ghost cell-center
             fluxx(i,j,k) = mu * ( phi(i,j,k) - phi(i-1,j,k) ) / (0.5d0*dx)
          end do
       end do
       !$omp end parallel do
    else if (adv_bc(1,1) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxx(lo(1),lo(2):hi(2),lo(3):hi(3)) = 0.d0
    end if

    ! hi-x boundary conditions
    if (adv_bc(1,2) .eq. EXT_DIR) then
       i=hi(1)+1
       !$omp parallel do private(j,k)
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             ! divide by 0.5*dx since the ghost cell value represents
             ! the value at the wall, not the ghost cell-center
             fluxx(i,j,k) = mu * ( phi(i,j,k) - phi(i-1,j,k) ) / (0.5d0*dx)
          end do
       end do
       !$omp end parallel do
    else if (adv_bc(1,2) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxx(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = 0.d0
    end if

    ! y-fluxes
    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             fluxy(i,j,k) = mu * ( phi(i,j,k) - phi(i,j-1,k) ) / dx
          end do
       end do
    end do
    !$omp end parallel do

    ! lo-y boundary conditions
    if (adv_bc(2,1) .eq. EXT_DIR) then
       j=lo(2)
       !$omp parallel do private(i,k)
       do k=lo(3),hi(3)
          do i=lo(1),hi(1)
             ! divide by 0.5*dx since the ghost cell value represents
             ! the value at the wall, not the ghost cell-center
             fluxy(i,j,k) = mu * ( phi(i,j,k) - phi(i,j-1,k) ) / (0.5d0*dx)
          end do
       end do
       !$omp end parallel do
    else if (adv_bc(2,1) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxy(lo(1):hi(1),lo(2),lo(3):hi(3)) = 0.d0
    end if

    ! hi-y boundary conditions
    if (adv_bc(2,2) .eq. EXT_DIR) then
       j=hi(2)+1
       !$omp parallel do private(i,k)
       do k=lo(3),hi(3)
          do i=lo(1),hi(1)
             ! divide by 0.5*dx since the ghost cell value represents
             ! the value at the wall, not the ghost cell-center
             fluxy(i,j,k) = mu * ( phi(i,j,k) - phi(i,j-1,k) ) / (0.5d0*dx)
          end do
       end do
       !$omp end parallel do
    else if (adv_bc(2,2) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxy(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = 0.d0
    end if

    ! z-fluxes
    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             fluxz(i,j,k) = mu * ( phi(i,j,k) - phi(i,j,k-1) ) / dx
          end do
       end do
    end do
    !$omp end parallel do

    ! lo-z boundary conditions
    if (adv_bc(3,1) .eq. EXT_DIR) then
       k=lo(3)
       !$omp parallel do private(i,j)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! divide by 0.5*dx since the ghost cell value represents
             ! the value at the wall, not the ghost cell-center
             fluxz(i,j,k) = mu * ( phi(i,j,k) - phi(i,j,k-1) ) / (0.5d0*dx)
          end do
       end do
       !$omp end parallel do
    else if (adv_bc(3,1) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxz(lo(1):hi(1),lo(2):lo(3),lo(3)) = 0.d0
    end if

    ! hi-z boundary conditions
    if (adv_bc(3,2) .eq. EXT_DIR) then
       k=hi(3)+1
       !$omp parallel do private(i,j)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! divide by 0.5*dx since the ghost cell value represents
             ! the value at the wall, not the ghost cell-center
             fluxz(i,j,k) = mu * ( phi(i,j,k) - phi(i,j,k-1) ) / (0.5d0*dx)
          end do
       end do
       !$omp end parallel do
    else if (adv_bc(3,2) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxz(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = 0.d0
    end if

    ! ****************************************************************
    ! Add the advective fluxes (here we assume constant velocity so no bc's)
    ! ****************************************************************

    ! x-fluxes: add the advective fluxes -- these are 2nd order slopes
    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1

             ! Assume that uadv > 0, which means we need the slope in cell (i-1,j,k)
             dlft = 2.d0*(phi(i-1,j,k) - phi(i-2,j,k))
             drgt = 2.d0*(phi(i  ,j,k) - phi(i-1,j,k))
             dcen = .25d0 * (dlft+drgt)
             if (dlft*drgt .ge. 0.d0) then
                dlim = min( abs(dlft), abs(drgt) )
             else
                dlim = 0.d0
             endif
             slope = sign(1.d0, dcen) * min( dlim, abs(dcen) )
          
             phi_edge = phi(i-1,j,k) + 0.5d0 * (1.d0 - uadv*dt/dx)*slope
   
             fluxx(i,j,k) = fluxx(i,j,k) - uadv * phi_edge

          end do
       end do
    end do
    !$omp end parallel do

    ! y-fluxes: add the advective fluxes
    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)

             ! Assume that vadv > 0, which means we need the slope in cell (i,j-1,k)
             dlft = 2.d0*(phi(i,j-1,k) - phi(i,j-2,k))
             drgt = 2.d0*(phi(i,j  ,k) - phi(i,j-1,k))
             dcen = .25d0 * (dlft+drgt)
             if (dlft*drgt .ge. 0.d0) then
                dlim = min( abs(dlft), abs(drgt) )
             else
                dlim = 0.d0
             endif
             slope = sign(1.d0, dcen) * min( dlim, abs(dcen) )
          
             phi_edge = phi(i,j-1,k) + 0.5d0 * (1.d0 - vadv*dt/dx) * slope
   
             fluxy(i,j,k) = fluxy(i,j,k) - vadv * phi_edge
   
          end do
       end do
    end do
    !$omp end parallel do

    ! z-fluxes: add the advective fluxes
    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! Assume that wadv > 0, which means we need the slope in cell (i,j,k-1)
             dlft = 2.d0*(phi(i,j,k-1) - phi(i,j,k-2))
             drgt = 2.d0*(phi(i,j,k  ) - phi(i,j,k-1))
             dcen = .25d0 * (dlft+drgt)
             if (dlft*drgt .ge. 0.d0) then
                dlim = min( abs(dlft), abs(drgt) )
             else
                dlim = 0.d0
             endif
             slope = sign(1.d0, dcen) * min( dlim, abs(dcen) )
             
             phi_edge = phi(i,j,k-1) + 0.5d0 * (1.d0 - wadv*dt/dx) * slope
   
             fluxz(i,j,k) = fluxz(i,j,k) - wadv * phi_edge
   
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine compute_flux_3d

end module compute_flux_module

