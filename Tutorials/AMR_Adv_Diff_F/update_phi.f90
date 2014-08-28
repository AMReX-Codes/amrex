
module update_phi_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use multifab_fill_ghost_module
  use ml_cc_restriction_module
  use bndry_reg_module

  implicit none

  private

  public :: update_phi, update_phi_single_level

contains

  subroutine update_phi(mla,phi_old,phi_new,flux,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi_old(:)
    type(multifab) , intent(inout) :: phi_new(:)
    type(multifab) , intent(in   ) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: dm, ng_p, ng_f, n, nlevs

    dm    = mla%dim
    nlevs = mla%nlevel

    ng_p = phi_old(1)%ng
    ng_f = flux(1,1)%ng

    do n=1,nlevs

       call update_phi_single_level(mla,phi_old(n),phi_new(n),flux(n,:),dx(n))

    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(phi_new(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(phi_new(nlevs),1,1,1,the_bc_tower%bc_tower_array(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1
          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(phi_new(n-1),phi_new(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(phi_new(n),phi_new(n-1),ng_p,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n), &
                                         1,1,1)
       end do

    end if

  end subroutine update_phi

  subroutine update_phi_single_level(mla,phi_old,phi_new,flux,dx)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi_old
    type(multifab) , intent(inout) :: phi_new
    type(multifab) , intent(in   ) :: flux(:)
    real(kind=dp_t), intent(in   ) :: dx

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, ng_f, i

    real(kind=dp_t), pointer ::  po(:,:,:,:)
    real(kind=dp_t), pointer ::  pn(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)

    dm    = mla%dim

    ng_p = phi_old%ng
    ng_f = flux(1)%ng

    do i=1,nfabs(phi_old)
       po  => dataptr(phi_old,i)
       pn  => dataptr(phi_new,i)
       fxp => dataptr(flux(1),i)
       fyp => dataptr(flux(2),i)
       lo = lwb(get_box(phi_old,i))
       hi = upb(get_box(phi_old,i))
       select case(dm)
       case (2)
          call update_phi_2d(po(:,:,1,1), pn(:,:,1,1), ng_p, &
                             fxp(:,:,1,1), fyp(:,:,1,1), ng_f, &
                             lo, hi, dx)
       case (3)
          fzp => dataptr(flux(3),i)
          call update_phi_3d(po(:,:,:,1), pn(:,:,:,1), ng_p, &
                             fxp(:,:,:,1), fyp(:,:,:,1), fzp(:,:,:,1), ng_f, &
                             lo, hi, dx)
       end select
    end do

  end subroutine update_phi_single_level

  subroutine update_phi_2d(phi_old, phi_new, ng_p, fluxx, fluxy, ng_f, lo, hi, dx)

    integer          :: lo(2), hi(2), ng_p, ng_f
    double precision :: phi_old(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision :: phi_new(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision ::   fluxx(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision ::   fluxy(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: dx

    ! local variables
    integer i,j

    ! Note that the factor of dt is already included in the fluxes
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          phi_new(i,j) = phi_old(i,j) + &
               ( fluxx(i+1,j)-fluxx(i,j) + fluxy(i,j+1)-fluxy(i,j) ) / dx
       end do
    end do

  end subroutine update_phi_2d

  subroutine update_phi_3d(phi_old, phi_new, ng_p, fluxx, fluxy, fluxz, ng_f, lo, hi, dx)

    integer          :: lo(3), hi(3), ng_p, ng_f
    double precision :: phi_old(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision :: phi_new(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision ::   fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision ::   fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision ::   fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: dx

    ! local variables
    integer i,j,k

    !$omp parallel do private(i,j,k)
    ! Note that the factor of dt is already included in the fluxes
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             phi_new(i,j,k) = phi_old(i,j,k) + &
                  ( fluxx(i+1,j,k)-fluxx(i,j,k) &
                   +fluxy(i,j+1,k)-fluxy(i,j,k) &
                   +fluxz(i,j,k+1)-fluxz(i,j,k) ) / dx
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine update_phi_3d

end module update_phi_module
