module init_data_module

  use multifab_module
  use layout_module

  implicit none

  private

  public :: init_data

contains
  
  subroutine init_data(level,data,dx,prob_lo)

    integer        , intent(in   ) :: level
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

    do i=1,nboxes(data)
       if ( multifab_remote(data,i) ) cycle
       dp => dataptr(data,i)
       lo = lwb(get_box(data,i))
       hi = upb(get_box(data,i))
       select case(dm)
       case (2)
          call init_data_2d(dp(:,:,1,1), ng, lo, hi, prob_lo, dx)
       case (3)
          call init_data_3d(dp(:,:,:,1), ng, lo, hi, prob_lo, dx)
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
    double precision :: U(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    double precision :: prob_lo(2)
    double precision :: dx
 
    ! local varables
    integer          :: i,j
    double precision :: x,y,dist

    do j = lo(2),hi(2)
       y = (dble(j)+0.5d0) * dx - 0.5d0

       do i = lo(1),hi(1)
          x = (dble(i)+0.5d0) * dx - 0.5d0

          dist = sqrt(x*x + y*y)
          U(i,j) = 0.5d0 * (1.d0 - tanh((dist-0.2d0)/0.025d0))

         end do
    end do

    end subroutine init_data_2d

    subroutine init_data_3d(U, ng, lo, hi, prob_lo, dx)

    integer          :: lo(3), hi(3), ng
    double precision :: U(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    double precision :: prob_lo(3)
    double precision :: dx
 
    ! local varables
    integer          :: i,j,k
    double precision :: x,y,z,dist

    do k = lo(3),hi(3)
       z = (dble(k)+0.5d0) * dx - 0.5d0

       do j = lo(2),hi(2)
          y = (dble(j)+0.5d0) * dx - 0.5d0

          do i = lo(1),hi(1)
             x = (dble(i)+0.5d0) * dx - 0.5d0

             dist = sqrt(x*x + y*y + z*z)
             U(i,j,k) = 0.5d0 * (1.d0 - tanh((dist-0.2d0)/0.025d0))

          end do
       end do
    end do

  end subroutine init_data_3d

end module init_data_module
