module init_phi_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use multifab_physbc_module
  use multifab_fill_ghost_module
  use ml_cc_restriction_module

  implicit none

  private

  public :: init_phi_on_level, init_phi

contains

  subroutine init_phi_on_level(phi,dx,prob_lo,the_bc_level)

    type(multifab) , intent(inout) :: phi
    real(kind=dp_t), intent(in   ) :: dx
    real(kind=dp_t), intent(in   ) :: prob_lo(phi%dim)
    type(bc_level) , intent(in   ) :: the_bc_level

    ! local
    integer i,ng,dm
    integer :: lo(phi%dim), hi(phi%dim)

    real(kind=dp_t), pointer :: dp(:,:,:,:)

    ng = phi%ng
    dm = phi%dim

    do i=1,nfabs(phi)
       dp => dataptr(phi,i)
       lo = lwb(get_box(phi,i))
       hi = upb(get_box(phi,i))
       select case(dm)
       case (2)
          call init_phi_2d(dp(:,:,1,1), ng, lo, hi, prob_lo, dx)
       case (3)
          call init_phi_3d(dp(:,:,:,1), ng, lo, hi, prob_lo, dx)
       end select
    end do

    call multifab_fill_boundary(phi)
    
    call multifab_physbc(phi,1,1,1,the_bc_level)

  end subroutine init_phi_on_level
  
  subroutine init_phi(mla,phi,dx,prob_lo,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: prob_lo(mla%dim)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: nlevs, dm, ng, i, n

    real(kind=dp_t), pointer :: dp(:,:,:,:)

    ng = phi(1)%ng
    dm = mla%dim
    nlevs = mla%nlevel

    do n=1,nlevs

       do i=1,nfabs(phi(n))
          dp => dataptr(phi(n),i)
          lo = lwb(get_box(phi(n),i))
          hi = upb(get_box(phi(n),i))
          select case(dm)
          case (2)
             call init_phi_2d(dp(:,:,1,1), ng, lo, hi, prob_lo, dx(n))
          case (3)
             call init_phi_3d(dp(:,:,:,1), ng, lo, hi, prob_lo, dx(n))
          end select
       end do
       
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
          call multifab_fill_ghost_cells(phi(n),phi(n-1),phi(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n), &
                                         1,1,1)
       end do

    end if

  end subroutine init_phi

  subroutine init_phi_2d(phi, ng, lo, hi, prob_lo, dx)

    integer          :: lo(2), hi(2), ng
    double precision :: phi(lo(1)-ng:,lo(2)-ng:)
    double precision :: prob_lo(2)
    double precision :: dx
 
    ! local varables
    integer          :: i,j
    double precision :: x,y,r1,r2,r3

    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+0.5d0) * dx
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+0.5d0) * dx

          r1 = ((x-0.5d0)**2 + (y-0.5d0)**2) / 0.01d0
          r2 = ((x-0.0d0)**2 + (y-0.0d0)**2) / 0.01d0
          r3 = ((x+0.5d0)**2 + (y+0.5d0)**2) / 0.01d0

          phi(i,j) = 1.d0 + exp(-r1) + exp(-r2) + exp(-r3)

       end do
    end do

    end subroutine init_phi_2d

    subroutine init_phi_3d(phi, ng, lo, hi, prob_lo, dx)

    integer          :: lo(3), hi(3), ng
    double precision :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    double precision :: prob_lo(3)
    double precision :: dx
 
    ! local varables
    integer          :: i,j,k
    double precision :: x,y,z,r1,r2,r3

    !$omp parallel do private(i,j,k,x,y,z,r2)
    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+0.5d0) * dx
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0) * dx
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0) * dx

             r1 = ((x-0.5d0)**2 + (y-0.5d0)**2 + (z-0.0d0)**2) / 0.01d0
             r2 = ((x-0.0d0)**2 + (y-0.0d0)**2 + (z-0.0d0)**2) / 0.01d0
             r3 = ((x+0.5d0)**2 + (y+0.5d0)**2 + (z-0.0d0)**2) / 0.01d0

             phi(i,j,k) = 1.d0 + exp(-r1) + exp(-r2) + exp(-r3)

          end do
       end do
    end do
    !$omp end parallel do

  end subroutine init_phi_3d

end module init_phi_module
