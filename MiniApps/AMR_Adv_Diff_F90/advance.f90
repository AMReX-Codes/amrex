
module advance_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use multifab_fill_ghost_module
  use ml_cc_restriction_module
  use bndry_reg_module
  use compute_flux_module
  use update_phi_module

  implicit none

  private

  public :: advance

contains

  subroutine advance(mla,phi_old,phi_new,bndry_flx,dx,dt,the_bc_tower,do_subcycling,num_substeps)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi_old(:)
    type(multifab) , intent(inout) :: phi_new(:)
    type(bndry_reg), intent(inout) :: bndry_flx(2:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    logical        , intent(in   ) :: do_subcycling
    integer        , intent(in   ) :: num_substeps

    ! local variables
    integer   :: i,istep,dm,n,nlevs,ng_p
    ! Array of edge-based multifabs; one for each direction
    type(multifab) :: flux(mla%dim,mla%nlevel)

    dm    = mla%dim
    nlevs = mla%nlevel
    ng_p = phi_new(1)%ng

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Explicit time advancement
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (do_subcycling) then

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! subcycling in time
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          n = 1
          istep = 1
          ! note that update_level is recursive
          call update_level(n,mla,phi_old,phi_new,bndry_flx,&
               dx,dt,the_bc_tower,istep,num_substeps)

    else

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! no subcycling in time
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! flux(i,n) has one component, zero ghost cells, and is nodal in direction i,
          !    i.e. it lives on the i-edges.
          do n = 1,nlevs
             do i = 1, dm
                call multifab_build_edge(flux(i,n),mla%la(n),1,0,i)
             end do
          end do
          
          ! Copy phi_new from the previous time step into phi_old for this time step
          do n = 1, nlevs
             call multifab_copy(mdst=phi_old(n),msrc=phi_new(n),ng=ng_p)
          end do

          ! Compute the face-centered gradients in each direction
          call compute_flux(mla,phi_old,flux,dx,dt,the_bc_tower)

          ! update phi using forward Euler discretization
          call update_phi(mla,phi_old,phi_new,flux,dx,dt,the_bc_tower)

          do n = 1,nlevs
             do i = 1, dm
                call multifab_destroy(flux(i,n))
             end do
          end do

    end if

  end subroutine advance

  recursive subroutine update_level(n,mla,phi_old,phi_new,bndry_flx,&
                                    dx,dt,the_bc_tower,step,num_substeps)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi_old(:)
    type(multifab) , intent(inout) :: phi_new(:)
    type(bndry_reg), intent(inout) :: bndry_flx(2:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: n,step
    integer        , intent(in   ) :: num_substeps

    real(kind=dp_t) :: alpha, scale
    integer         :: istep, i, dm, ng_p
    ! Array of edge-based multifabs; one for each direction
    type(multifab) :: flux(mla%dim)

    dm = mla%dim
    ng_p = phi_new(n)%ng

    ! We only want one processor to write to screen
    if ( parallel_IOProcessor() ) &
       print*,'   Advancing level: ',n,' with dt = ',dt(n)

    ! Copy phi_new from the previous time step into phi_old for this time step
    call multifab_copy(mdst=phi_old(n),msrc=phi_new(n),ng=ng_p)

    if (n .eq. 1) then 

        ! Fill ghost cells for two adjacent grids at the same level
        ! this includes periodic domain boundary ghost cells
        call multifab_fill_boundary(phi_old(n))
 
        ! Fill non-periodic domain boundary ghost cells
        call multifab_physbc(phi_old(n),1,1,1,the_bc_tower%bc_tower_array(n))

     else if (step .eq. 1) then
        ! First update boundary conditions
        ! Fill level n ghost cells using interpolation from level n-1 data
        ! Note that multifab_fill_boundary and multifab_physbc are called for
        ! both levels n-1 and n
        ! The coarsest level
        call multifab_fill_ghost_cells(phi_old(n),phi_old(n-1),ng_p,&
                                       mla%mba%rr(n-1,:), &
                                       the_bc_tower%bc_tower_array(n-1), &
                                       the_bc_tower%bc_tower_array(n), &
                                       1,1,1)
    else
        !
        ! This version of multifab_fill_ghost_cells takes both a crse_old and a crse_new
        !     and does interpolation in time as well as in space.  The coefficient alpha
        !     is used to define:   fine = interp (alpha*crse_old + (1-alpha)*crse_new )
        !      
        ! If num_substeps = 2, then at step = 1, alpha = 1.0
        !                           at step = 2, alpha = 0.5
        ! If num_substeps = 4, then at step = 1, alpha = 1.0
        !                           at step = 2, alpha = 0.75
        !                           at step = 3, alpha = 0.5
        !                           at step = 4, alpha = 0.25
        !      
        alpha = 1.d0 - dble(step-1)/dble(num_substeps)
        call multifab_fill_ghost_cells_t(phi_old(n),phi_old(n-1),phi_new(n-1),alpha,&
                                         ng_p,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n), &
                                         1,1,1)
    end if

    do i = 1, dm
       call multifab_build_edge(flux(i),mla%la(n),1,0,i)
    end do

    ! Compute fluxes at level n.
    call compute_flux_single_level(mla,phi_old(n),flux,dx(n),dt(n),the_bc_tower%bc_tower_array(n))

    ! Update solution at level n.
    call update_phi_single_level(mla,phi_old(n),phi_new(n),flux,dx(n),dt(n))

    if (n .gt. 1) then
       ! Copy fine fluxes from cell boundaries into boundary registers for use in refluxing.
       scale = dt(n)*dx(n)**(dm-1) ! dt*area
       call flux_reg_fine_add(bndry_flx(n),flux,scale)
    end if

    if (n .eq. mla%nlevel) then ! finest level

       do i = 1, dm
          call multifab_destroy(flux(i))
       end do

    else  ! Now recursively update the next finer level

       scale = -dt(n)*dx(n)**(dm-1) ! -dt*area
       call flux_reg_crse_init(bndry_flx(n+1), flux, scale)

       do i = 1, dm
          call multifab_destroy(flux(i))
       end do

       do istep = 1, num_substeps
          call update_level(n+1,mla,phi_old,phi_new,bndry_flx,&
               dx,dt,the_bc_tower,istep,num_substeps)
       end do

       if ( parallel_IOProcessor() ) &
            print*,'   Refluxing level: ',n

       ! note that bndry_flx contains sum_{fine}{flux*dt*area} - {flux*dt*area}_crse 
       scale = 1.d0/dx(n)**dm ! 1/volume 
       call reflux(phi_new(n),bndry_flx(n+1),scale)

       ! Average the level (n+1) phi down onto level n. This overwrites all
       !     previous values of phi at level n under level (n+1) grids.
       call ml_cc_restriction_c(phi_new(n),1,phi_new(n+1),1,mla%mba%rr(n,:),1)
    end if

  end subroutine update_level

end module advance_module

