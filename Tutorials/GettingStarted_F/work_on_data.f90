module work_on_data_module

  use multifab_module

  implicit none

  private

  public :: work_on_data

contains
  
  subroutine work_on_data(data)

    type(multifab), intent(inout) :: data

    integer                  :: i,dm,ng,nc
    real(kind=dp_t), pointer :: dp(:,:,:,:)

    integer, allocatable     :: lo(:),hi(:)

    dm = data%dim   ! dm is dimensionality
    ng = data%ng    ! ng is number of ghost cells
    nc = data%nc    ! nc is number of components

    allocate(lo(dm),hi(dm))

    ! loop over the grids owned by this processor
    do i=1,nfabs(data)
       dp => dataptr(data,i)          ! dp points to data inside fab
       lo = lwb(get_box(data,i))      ! get lo indices of box
       hi = upb(get_box(data,i))      ! get hi indices of box
       select case(dm)
       case (2)
          call work_on_data_2d(dp(:,:,1,:), ng, nc, lo, hi)
       case (3)
          call work_on_data_3d(dp(:,:,:,:), ng, nc, lo, hi)
       end select
    end do
    
    ! fill periodic domain boundary and neighboring grid ghost cells
    call multifab_fill_boundary(data)  

    deallocate(lo,hi)
    
  end subroutine work_on_data
  
  subroutine work_on_data_2d(data, ng, nc, lo, hi)

    integer          :: lo(2), hi(2), ng, nc
    double precision :: data(lo(1)-ng:,lo(2)-ng:,:)

    ! local variables
    integer :: i,j,n

    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
    do n=1,nc
       ! some silly function I made up
       data(i,j,n) = (i + j) * n
    end do
    end do
    end do

  end subroutine work_on_data_2d
  
  subroutine work_on_data_3d(data, ng, nc, lo, hi)

    integer          :: lo(3), hi(3), ng, nc
    double precision :: data(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)

    ! local variables
    integer :: i,j,k,n

    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
    do n=1,nc
       ! some silly function I made up
       data(i,j,k,n) = (i + j) * n
    end do
    end do
    end do
    end do

  end subroutine work_on_data_3d

end module work_on_data_module

