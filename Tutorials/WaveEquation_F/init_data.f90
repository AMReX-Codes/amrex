module init_data_module

  implicit none

  private

  public :: init_data

contains
  
  subroutine init_data(data,dx,prob_lo)

    use multifab_module

    type(multifab) , intent(inout) :: data
    real(kind=dp_t), intent(in   ) :: dx
    real(kind=dp_t), intent(in   ) :: prob_lo(data%dim)

    ! local variables
    integer :: lo(data%dim), hi(data%dim)
    integer :: dm, ng, i

    real(kind=dp_t), pointer :: dp(:,:,:,:)

    ! set these here so we don't have to pass them into the subroutine
    dm = data%dim
    ng = data%ng

    do i=1,nfabs(data)
       dp => dataptr(data,i)
       lo = lwb(get_box(data,i))
       hi = upb(get_box(data,i))
       select case(dm)
       case (2)
          call init_data_2d(dp(:,:,1,:), ng, lo, hi, prob_lo, dx)
       case (3)
          call init_data_3d(dp(:,:,:,:), ng, lo, hi, prob_lo, dx)
       end select
    end do

    ! fill ghost cells
    ! this only fills periodic ghost cells and ghost cells for neighboring
    ! grids at the same level.  Physical boundary ghost cells are filled
    ! using multifab_physbc.  But this problem is periodic, so this
    ! call is sufficient.
    call multifab_fill_boundary(data)

  end subroutine init_data

  subroutine init_data_2d(U, ng, lo, hi, prob_lo, dx)

    integer          :: lo(2), hi(2), ng
    double precision :: U(lo(1)-ng:,lo(2)-ng:,:)
    double precision :: prob_lo(2)
    double precision :: dx
 
    ! local varables
    integer          :: i,j
    double precision :: x,y,r2

    do j = lo(2), hi(2)
         y = prob_lo(2) + (dble(j)+0.5d0) * dx
         do i = lo(1), hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx
            r2 = (x*x + y*y) / 0.01

            U(i,j,1) = 0.d0
            U(i,j,2) = exp(-r2)

         end do
      end do

    end subroutine init_data_2d

    subroutine init_data_3d(U, ng, lo, hi, prob_lo, dx)

    integer          :: lo(3), hi(3), ng
    double precision :: U(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    double precision :: prob_lo(3)
    double precision :: dx
 
    ! local varables
    integer          :: i,j,k
    double precision :: x,y,z,r2

    do k = lo(3), hi(3)
       z = prob_lo(3) + (dble(k)+0.5d0) * dx
       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0) * dx
          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0) * dx
             r2 = (x*x + y*y + z*z) / 0.01

             U(i,j,k,1) = 0.d0
             U(i,j,k,2) = exp(-r2)

          end do
       end do
    end do

  end subroutine init_data_3d

end module init_data_module
