module init_data_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: init_data

contains
  
  subroutine init_data(mla,dx,prob_lo,data)

    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), intent(in   ) :: prob_lo(mla%dim)
    real(kind=dp_t), intent(in   ) :: dx(mla%nlevel,mla%dim)
    type(multifab) , intent(inout) :: data(mla%nlevel)

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: nlevs, dm, ng, n, i

    real(kind=dp_t), pointer :: dp(:,:,:,:)

    ! set these here so we don't have to pass them into the subroutine
    nlevs = mla%nlevel
    dm    = mla%dim
    ng    = nghost(data(1))

    do n=1,nlevs

       do i=1,nboxes(data(n))
          if ( multifab_remote(data(n),i) ) cycle
          dp => dataptr(data(n),i)
          lo = lwb(get_box(data(n),i))
          hi = upb(get_box(data(n),i))
          select case(dm)
          case (2)
             call init_data_2d(dp(:,:,1,:), ng, lo, hi, prob_lo, dx(n,:))
          case (3)
             call init_data_3d(dp(:,:,:,:), ng, lo, hi, prob_lo, dx(n,:))
          end select
       end do

       ! fill ghost cells
       ! this only fills periodic ghost cells and ghost cells for neighboring
       ! grids at the same level.  Physical boundary ghost cells are filled
       ! using multifab_physbc.  But this problem is periodic, so this
       ! call is sufficient.
       call multifab_fill_boundary(data(n))

    end do

  end subroutine init_data

  subroutine init_data_2d(U, ng, lo, hi, prob_lo, dx)

    integer          :: lo(2), hi(2), ng
    double precision :: U(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,2)
    double precision :: prob_lo(2)
    double precision :: dx(2)
 
    ! local varables
    integer          :: i,j
    double precision :: x,y,r2

    do j = lo(2), hi(2)
         y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
         do i = lo(1), hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
            r2 = (x*x + y*y) / 0.01

            U(i,j,1) = 0.d0
            U(i,j,2) = exp(-r2)

         end do
      end do

    end subroutine init_data_2d

    subroutine init_data_3d(U, ng, lo, hi, prob_lo, dx)

    integer          :: lo(3), hi(3), ng
    double precision :: U(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,2)
    double precision :: prob_lo(3)
    double precision :: dx(3)
 
    ! local varables
    integer          :: i,j,k
    double precision :: x,y,z,r2

    do k = lo(3), hi(3)
       z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
             r2 = (x*x + y*y + z*z) / 0.01

             U(i,j,k,1) = 0.d0
             U(i,j,k,2) = exp(-r2)

          end do
       end do
    end do

  end subroutine init_data_3d

end module init_data_module
