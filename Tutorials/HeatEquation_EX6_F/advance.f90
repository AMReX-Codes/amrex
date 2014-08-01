module advance_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use multifab_fill_ghost_module
  use ml_restriction_module
  use subcycling_module
  use bndry_reg_module

  implicit none

  private

  public :: advance

contains

  subroutine advance(mla,phi,isrefined,dx,dt,the_bc_tower,do_subcycling,num_substeps)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi(:)
    type(lmultifab), intent(in   ) :: isrefined(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    logical        , intent(in   ) :: do_subcycling
    integer        , intent(in   ) :: num_substeps

    ! local variables
    integer i,dm,n,nlevs,ng_p

    ! an array of multifabs; one for each direction
    type(multifab) :: flux(mla%nlevel,mla%dim)

    ! an array of multifabs to store phi from previous
    ! time step
    type(multifab) :: phi_prev_ts(mla%nlevel)
    ! an array of multifabs to store some temporary
    ! data needed while interpolating in time, only
    ! when (num_substeps.eq.4)
    type(multifab) :: delta_phi_in_time(mla%nlevel)
    !subcycling counters
    integer        :: lvl_computed_count(mla%nlevel)
    real(kind=dp_t):: dt_local

    type(multifab) :: alpha(mla%nlevel)
    type(multifab) :: beta(mla%nlevel,mla%dim)
    type(multifab) :: rhs(mla%nlevel)

    type(bndry_reg) :: fine_flx(2:mla%nlevel)

    real(dp_t) :: dx_vector(mla%nlevel,mla%dim)

    integer :: bc_comp

    dm = mla%dim
    nlevs = mla%nlevel

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! explicit time advancement
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (do_subcycling) then

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! subcycling in time
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ng_p = phi(1)%ng
          ! build the flux(:,:) multifabs
          do n=1,nlevs
            if (n.ne.1) call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n),1)
            do i=1,dm
                ! flux(n,i) has one component, zero ghost cells, and is nodal in direction i
                call multifab_build_edge(flux(n,i),mla%la(n),1,0,i)
            end do
          end do

          do n=1,nlevs
            ! build and initialize multifab storing phi from previous
            ! time step.
            ! Old values from cells overlayed by ghost cells
            ! are used to interpolate in time the state of ghost cells.
            call multifab_build(phi_prev_ts(n),mla%la(n),1,1)
            call multifab_build(delta_phi_in_time(n),mla%la(n),1,1)
            call setval(phi_prev_ts(n),0.0d0,all=.TRUE.)
            call setval(delta_phi_in_time(n),0.0d0,all=.TRUE.)
          end do

          ! initialize
          lvl_computed_count(:)=0
          ! start at coarse level (1)
          n=1
          ! we only want one processor to write to screen
          if ( parallel_IOProcessor() ) then
             print*,' advancing level: ',n
          end if

          ! fill ghost cells for two adjacent grids at the same level
          ! this includes periodic domain boundary ghost cells
          call multifab_fill_boundary(phi(n))
          ! fill non-periodic domain boundary ghost cells
          call multifab_physbc(phi(n),1,1,1,the_bc_tower%bc_tower_array(n))

          call compute_flux_single_level(n,mla,phi,flux,dx,the_bc_tower)
          ! store copy of data at this level for ghost cells
          ! TODO: do this only for cells overlayed by ghost cells)
          ! phi_prev_ts(n) = phi(n)
          if (nlevs .ne. 1) call multifab_copy(mdst=phi_prev_ts(n),msrc=phi(n),ng=ng_p)
          ! update solution at coarse level (n=1)
          dt_local = dt / (dble(num_substeps))**(n-1)
          call update_phi_single_level(n,mla,phi,flux,dx,dt_local)

          ! if more than one level, use AMR
          if (nlevs .ne. 1) then

             n=2
             ! start loop of time stepping procedure
             do while (n.ne.1)
               if (lvl_computed_count(n).lt.num_substeps) then
                 ! we only want one processor to write to screen
                 if ( parallel_IOProcessor() ) then
                    print*,' advancing level: ',n
                 end if
                 ! if n is not the highest level, store states for ghost cells of level n+1
                 ! this is the only way to access this data later,
                 ! since lower level is updated beforehand
                 if (n.ne.nlevs) call multifab_copy(mdst=phi_prev_ts(n),msrc=phi(n),ng=ng_p)
                 ! first update boundary conditions
                 ! fill level n ghost cells using interpolation from level n-1 data
                 ! note that multifab_fill_boundary and multifab_physbc are called for
                 ! both levels n-1 and n
                 if (lvl_computed_count(n).eq.0) then
                   ! this case applies when (num_substeps.eq.2) and (num_substeps.eq.4)
                   ! apply ghost cells state from previous coarse time (sub-)step
                   call multifab_fill_ghost_cells(phi(n),phi_prev_ts(n-1),ng_p,&
                         mla%mba%rr(n-1,:), &
                         the_bc_tower%bc_tower_array(n-1), &
                         the_bc_tower%bc_tower_array(n), &
                         1,1,1)
                 else if ((lvl_computed_count(n).eq.1).and.(num_substeps.eq.2)) then
                   ! this case only applies when (num_substeps.eq.2)
                   ! we use simplified algorithm, without use of delta_phi_in_time multifab
                   ! (cf. 'else if' below)
                   ! phi_prev_ts(n-1) = 0.5*(phi(n-1)+phi_prev_ts(n-1))
                   call multifab_copy(mdst=phi_prev_ts(n-1),msrc=phi(n-1),ng=ng_p,&
                                      filter=filter_mean)
                   ! apply ghost cells state interpolated in time
                   call multifab_fill_ghost_cells(phi(n),phi_prev_ts(n-1),ng_p,&
                         mla%mba%rr(n-1,:), &
                         the_bc_tower%bc_tower_array(n-1), &
                         the_bc_tower%bc_tower_array(n), &
                         1,1,1)
                 else if ((lvl_computed_count(n).eq.1).and.(num_substeps.eq.4)) then
                   ! this case only applies when (num_substeps.eq.4)
                   ! delta_phi_in_time(n-1) = 0.25*(phi_prev_ts(n-1)-phi(n-1))
                   call multifab_copy(mdst=delta_phi_in_time(n-1),msrc=phi_prev_ts(n-1),ng=ng_p,&
                                      filter=filter_assign)
                   call multifab_copy(mdst=delta_phi_in_time(n-1),msrc=phi(n-1),ng=ng_p,&
                                      filter=filter_quarter_of_diff)
                   ! phi_prev_ts(n-1) = phi_prev_ts(n-1)) + delta_phi_in_time(n-1)
                   call multifab_copy(mdst=phi_prev_ts(n-1),msrc=delta_phi_in_time(n-1),ng=ng_p,&
                                      filter=filter_add)
                   ! apply ghost cells state interpolated in time
                   call multifab_fill_ghost_cells(phi(n),phi_prev_ts(n-1),ng_p,&
                         mla%mba%rr(n-1,:), &
                         the_bc_tower%bc_tower_array(n-1), &
                         the_bc_tower%bc_tower_array(n), &
                         1,1,1)
                 else ! if ((lvl_computed_count(n).ge.2).and.(num_substeps.eq.4)) then
                   ! this case only applies when (num_substeps.eq.4) and when
                   ! (lvl_computed_count(n).eq.2) or (lvl_computed_count(n).eq.3)
                   ! phi_prev_ts(n-1) = phi_prev_ts(n-1)) + delta_phi_in_time(n-1)
                   call multifab_copy(mdst=phi_prev_ts(n-1),msrc=delta_phi_in_time(n-1),ng=ng_p,&
                                      filter=filter_add)
                   ! apply ghost cells state interpolated in time
                   call multifab_fill_ghost_cells(phi(n),phi_prev_ts(n-1),ng_p,&
                         mla%mba%rr(n-1,:), &
                         the_bc_tower%bc_tower_array(n-1), &
                         the_bc_tower%bc_tower_array(n), &
                         1,1,1)
                 endif
                 ! and then compute fluxes at level n
                 call compute_flux_single_level(n,mla,phi,flux,dx,the_bc_tower)
                 ! update solution at level n
                 dt_local = dt / (dble(num_substeps))**(n-1)
                 call update_phi_single_level(n,mla,phi,flux,dx,dt_local)
                 ! store fine fluxes for flux fix

                 ! set sum of fine fluxes to 0 at the first substep at this level
                 if (lvl_computed_count(n).eq.0) call bndry_reg_setval(fine_flx(n), 0.0d0)

                 ! store fine fluxes from grid boundaries to boundary registers
                 call flux_bndry_reg_copy_c(fine_flx(n),1,flux(n,:),1,1,filter=filter_add)

                 ! if this is the last substep, restrict fine fluxes from grid boundaries and substract
                 ! these from the flux at the lower level
                 ! coarse_flux = coarse_flux - 1/(num_substeps) * sum_of_restricted_fine_fluxes
                 if (lvl_computed_count(n).eq.(num_substeps-1)) then
                    if (num_substeps.eq.2) then
                       call fine_flx_restriction_c(flux(n-1,:),1,fine_flx(n),1,mla%mba%rr(n-1,:),1,&
                                                   filter=filter_subtract_half)
                    else !if (num_substeps.eq.4) then
                       call fine_flx_restriction_c(flux(n-1,:),1,fine_flx(n),1,mla%mba%rr(n-1,:),1,&
                                                   filter=filter_subtract_quarter)
                    endif
                 endif

                 ! count
                 lvl_computed_count(n)=lvl_computed_count(n)+1
                 ! if n is the highest level, stay on this level, othewise go to higher one
                 if (n.eq.nlevs) then
                   !n=n
                 else
                   n=n+1
                 endif
               else
                 ! reset counter
                 lvl_computed_count(n)=0
                 ! restriction (downsampling)
                 call ml_cc_restriction_c(phi(n-1),1,phi(n),1,mla%mba%rr(n-1,:),1)
                 ! go to previous level
                 n=n-1
                 ! update phi by applying flux fix
                 dt_local = dt / (dble(num_substeps))**(n-1)
                 call apply_fluxfix_single_level(n,mla,phi,flux,isrefined,dx,dt_local)
               endif
             end do

          endif

          ! make sure to destroy the multifabs or you'll leak memory
          do n=1,nlevs
            if (n.ne.1) call bndry_reg_destroy(fine_flx(n))
            do i=1,dm
                call multifab_destroy(flux(n,i))
            end do
            call multifab_destroy(phi_prev_ts(n))
            call multifab_destroy(delta_phi_in_time(n))
          end do

    else

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! no subcycling in time
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

    end if

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
    integer :: dm, ng_p, ng_f, n, nlevs

    dm    = mla%dim
    nlevs = mla%nlevel

    ng_p = phi(1)%ng
    ng_f = flux(1,1)%ng

    do n=1,nlevs

       call update_phi_single_level(n,mla,phi,flux,dx,dt)

    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(phi(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(phi(nlevs),1,1,1,the_bc_tower%bc_tower_array(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1
          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(phi(n-1),phi(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(phi(n),phi(n-1),ng_p,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n), &
                                         1,1,1)
       end do

    end if

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

