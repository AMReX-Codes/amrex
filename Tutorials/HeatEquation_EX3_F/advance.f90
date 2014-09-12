module advance_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use ml_cc_restriction_module
  use ml_restrict_fill_module

  implicit none

  private

  public :: advance

contains
  
  subroutine advance(mla,phi,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer i,dm,n,nlevs

    ! an array of multifabs; one for each direction
    type(multifab) :: flux(mla%nlevel,mla%dim) 

    dm = mla%dim
    nlevs = mla%nlevel

    ! build the flux(:,:) multifabs
    do n=1,nlevs
       do i=1,dm
          ! flux(n,i) has one component, zero ghost cells, and is nodal in direction i
          call multifab_build_edge(flux(n,i),mla%la(n),1,0,i)
       end do
    end do

    ! compute the face-centered gradients in each direction
    call compute_flux(mla,phi,flux,dx,the_bc_tower)
    
    ! update phi using forward Euler discretization
    call update_phi(mla,phi,flux,dx,dt,the_bc_tower)

    ! make sure to destroy the multifab or you'll leak memory
    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(flux(n,i))
       end do
    end do

  end subroutine advance

  subroutine compute_flux(mla,phi,flux,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi(:)
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: dm, ng_p, ng_f, i, n, nlevs

    dm    = mla%dim
    nlevs = mla%nlevel

    ng_p = phi(1)%ng
    ng_f = flux(1,1)%ng

    do n=1,nlevs

       call compute_flux_single_level(n,mla,phi,flux,dx,the_bc_tower)

    end do

    ! set level n-1 fluxes to be the average of the level n fluxes covering it
    ! the loop over nlevs must count backwards to make sure the finer grids are done first
    do n=nlevs,2,-1
       do i=1,dm
          call ml_edge_restriction_c(flux(n-1,i),1,flux(n,i),1,mla%mba%rr(n-1,:),i,1)
       end do
    end do

  end subroutine compute_flux

  subroutine compute_flux_single_level(n,mla,phi,flux,dx,the_bc_tower)

    integer        , intent(in   ) :: n
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi(:)
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, ng_f, i

    real(kind=dp_t), pointer ::  pp(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)

    dm    = mla%dim

    ng_p = phi(1)%ng
    ng_f = flux(1,1)%ng

    do i=1,nfabs(phi(n))
       pp  => dataptr(phi(n),i)
       fxp => dataptr(flux(n,1),i)
       fyp => dataptr(flux(n,2),i)
       lo = lwb(get_box(phi(n),i))
       hi = upb(get_box(phi(n),i))
       select case(dm)
       case (2)
          call compute_flux_2d(pp(:,:,1,1), ng_p, &
                               fxp(:,:,1,1),  fyp(:,:,1,1), ng_f, &
                               lo, hi, dx(n), &
                               the_bc_tower%bc_tower_array(n)%adv_bc_level_array(i,:,:,1))
       case (3)
          fzp => dataptr(flux(n,3),i)
          call compute_flux_3d(pp(:,:,:,1), ng_p, &
                               fxp(:,:,:,1),  fyp(:,:,:,1), fzp(:,:,:,1), ng_f, &
                               lo, hi, dx(n), &
                               the_bc_tower%bc_tower_array(n)%adv_bc_level_array(i,:,:,1))
       end select
    end do

  end subroutine compute_flux_single_level

  subroutine compute_flux_2d(phi, ng_p, fluxx, fluxy, ng_f, lo, hi, dx, adv_bc)

    integer          :: lo(2), hi(2), ng_p, ng_f
    double precision ::   phi(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: dx
    integer          :: adv_bc(:,:)

    ! local variables
    integer i,j

    ! x-fluxes
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / dx
       end do
    end do

    ! lo-x boundary conditions
    if (adv_bc(1,1) .eq. EXT_DIR) then
       i=lo(1)
       do j=lo(2),hi(2)
          ! divide by 0.5*dx since the ghost cell value represents
          ! the value at the wall, not the ghost cell-center
          fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / (0.5d0*dx)
       end do
    else if (adv_bc(1,1) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxx(lo(1),lo(2):hi(2)) = 0.d0
    end if

    ! hi-x boundary conditions
    if (adv_bc(1,2) .eq. EXT_DIR) then
       i=hi(1)+1
       do j=lo(2),hi(2)
          ! divide by 0.5*dx since the ghost cell value represents
          ! the value at the wall, not the ghost cell-center
          fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / (0.5d0*dx)
       end do
    else if (adv_bc(1,2) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxx(hi(1)+1,lo(2):hi(2)) = 0.d0
    end if

    ! y-fluxes
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / dx
       end do
    end do

    ! lo-y boundary conditions
    if (adv_bc(2,1) .eq. EXT_DIR) then
       j=lo(2)
       do i=lo(1),hi(1)
          ! divide by 0.5*dx since the ghost cell value represents
          ! the value at the wall, not the ghost cell-center
          fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / (0.5d0*dx)
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
          fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / (0.5d0*dx)
       end do
    else if (adv_bc(2,2) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxy(lo(1):hi(1),hi(2)+1) = 0.d0
    end if
    
  end subroutine compute_flux_2d

  subroutine compute_flux_3d(phi, ng_p, fluxx, fluxy, fluxz, ng_f, &
                             lo, hi, dx, adv_bc)

    integer          :: lo(3), hi(3), ng_p, ng_f
    double precision ::   phi(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: dx
    integer          :: adv_bc(:,:)

    ! local variables
    integer i,j,k

    ! x-fluxes
    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             fluxx(i,j,k) = ( phi(i,j,k) - phi(i-1,j,k) ) / dx
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
             fluxx(i,j,k) = ( phi(i,j,k) - phi(i-1,j,k) ) / (0.5d0*dx)
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
             fluxx(i,j,k) = ( phi(i,j,k) - phi(i-1,j,k) ) / (0.5d0*dx)
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
             fluxy(i,j,k) = ( phi(i,j,k) - phi(i,j-1,k) ) / dx
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
             fluxy(i,j,k) = ( phi(i,j,k) - phi(i,j-1,k) ) / (0.5d0*dx)
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
             fluxy(i,j,k) = ( phi(i,j,k) - phi(i,j-1,k) ) / (0.5d0*dx)
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
             fluxz(i,j,k) = ( phi(i,j,k) - phi(i,j,k-1) ) / dx
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
             fluxz(i,j,k) = ( phi(i,j,k) - phi(i,j,k-1) ) / (0.5d0*dx)
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
             fluxz(i,j,k) = ( phi(i,j,k) - phi(i,j,k-1) ) / (0.5d0*dx)
          end do
       end do
       !$omp end parallel do
    else if (adv_bc(3,2) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxz(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = 0.d0
    end if

  end subroutine compute_flux_3d

  subroutine update_phi(mla,phi,flux,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi(:)
    type(multifab) , intent(in   ) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: n, nlevs

    nlevs = mla%nlevel

    do n=1,nlevs

       call update_phi_single_level(n,mla,phi,flux,dx,dt)

    end do

    ! restrict the multi-level data, and
    ! fill all boundaries: same-level, coarse-fine, periodic, and domain boundaries 
    call ml_restrict_and_fill(nlevs, phi, mla%mba%rr, the_bc_tower%bc_tower_array)

  end subroutine update_phi

  subroutine update_phi_single_level(n,mla,phi,flux,dx,dt)

    integer        , intent(in   ) :: n
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi(:)
    type(multifab) , intent(in   ) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, ng_f, i

    real(kind=dp_t), pointer ::  pp(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)

    dm    = mla%dim

    ng_p = phi(1)%ng
    ng_f = flux(1,1)%ng

    do i=1,nfabs(phi(n))
       pp  => dataptr(phi(n),i)
       fxp => dataptr(flux(n,1),i)
       fyp => dataptr(flux(n,2),i)
       lo = lwb(get_box(phi(n),i))
       hi = upb(get_box(phi(n),i))
       select case(dm)
       case (2)
          call update_phi_2d(pp(:,:,1,1), ng_p, &
                             fxp(:,:,1,1), fyp(:,:,1,1), ng_f, &
                             lo, hi, dx(n), dt)
       case (3)
          fzp => dataptr(flux(n,3),i)
          call update_phi_3d(pp(:,:,:,1), ng_p, &
                             fxp(:,:,:,1), fyp(:,:,:,1), fzp(:,:,:,1), ng_f, &
                             lo, hi, dx(n), dt)
       end select
    end do

  end subroutine update_phi_single_level

  subroutine update_phi_2d(phi, ng_p, fluxx, fluxy, ng_f, lo, hi, dx, dt)

    integer          :: lo(2), hi(2), ng_p, ng_f
    double precision ::   phi(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: dx, dt

    ! local variables
    integer i,j

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          phi(i,j) = phi(i,j) + dt * &
               ( fluxx(i+1,j)-fluxx(i,j) + fluxy(i,j+1)-fluxy(i,j) ) / dx
       end do
    end do

  end subroutine update_phi_2d

  subroutine update_phi_3d(phi, ng_p, fluxx, fluxy, fluxz, ng_f, lo, hi, dx, dt)

    integer          :: lo(3), hi(3), ng_p, ng_f
    double precision ::   phi(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: dx, dt

    ! local variables
    integer i,j,k

    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             phi(i,j,k) = phi(i,j,k) + dt * &
                  ( fluxx(i+1,j,k)-fluxx(i,j,k) &
                   +fluxy(i,j+1,k)-fluxy(i,j,k) &
                   +fluxz(i,j,k+1)-fluxz(i,j,k) ) / dx
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine update_phi_3d

end module advance_module

