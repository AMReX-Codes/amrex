module init_phi_module

  use multifab_module
  use define_bc_module
  use multifab_physbc_module

  implicit none

  private

  public :: init_phi

contains
  
  subroutine init_phi(phi,dx,prob_lo,the_bc_tower)

    type(multifab) , intent(inout) :: phi
    real(kind=dp_t), intent(in   ) :: dx
    real(kind=dp_t), intent(in   ) :: prob_lo(phi%dim)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: lo(phi%dim), hi(phi%dim)
    integer :: dm, ng, i

    real(kind=dp_t), pointer :: dp(:,:,:,:)

    dm = phi%dim
    ng = phi%ng

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

    ! fill ghost cells for two adjacent grids at the same level
    ! this includes periodic domain boundary ghost cells
    call multifab_fill_boundary(phi)

    ! fill non-periodic domain boundary ghost cells
    call multifab_physbc(phi,1,1,1,the_bc_tower%bc_tower_array(1))

  end subroutine init_phi

  subroutine init_phi_2d(phi, ng, lo, hi, prob_lo, dx)

    integer          :: lo(2), hi(2), ng
    double precision :: phi(lo(1)-ng:,lo(2)-ng:)
    double precision :: prob_lo(2)
    double precision :: dx
 
    ! local varables
    integer          :: i,j
    double precision :: x,y,r2

    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+0.5d0) * dx
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+0.5d0) * dx

          r2 = ((x-0.25d0)**2 + (y-0.25d0)**2) / 0.01d0
          phi(i,j) = 1.d0 + exp(-r2)

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
    double precision :: x,y,z,r2

    !$omp parallel do private(i,j,k,x,y,z,r2)
    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+0.5d0) * dx
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0) * dx
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0) * dx

             r2 = ((x-0.25d0)**2 + (y-0.25d0)**2 + (z-0.25d0)**2) / 0.01d0
             phi(i,j,k) = 1.d0 + exp(-r2)

          end do
       end do
    end do
    !$omp end parallel do

  end subroutine init_phi_3d

end module init_phi_module
