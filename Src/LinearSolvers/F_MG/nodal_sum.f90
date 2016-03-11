module nodal_sum_module

  use multifab_module

  implicit none

  private

  public :: nodal_sum

contains

  ! returns the nodal sum and volume of level 1 only
  ! interior nodes get a weighting of "1"
  ! face nodes get a weighting of "1/2"
  ! corner (in 2D) and edge (in 3D) nodes get a weighting of "1/4"
  ! corner (in 3D) nodes get a weighting of "1/8"
  subroutine nodal_sum(mf,sum,vol)

    type(multifab) , intent(in   ) :: mf(:)
    real(kind=dp_t), intent(inout) :: sum, vol

    integer :: dm,i,ng

    real(kind=dp_t), pointer :: dp(:,:,:,:)
    integer :: lo(multifab_get_dim(mf(1)))
    integer :: hi(multifab_get_dim(mf(1)))

    real(kind=dp_t) :: sum_proc, vol_proc

    dm = multifab_get_dim(mf(1))
    ng = mf(1)%ng

    sum_proc = 0.d0
    vol_proc = 0.d0

    do i=1,nfabs(mf(1))
       dp => dataptr(mf(1),i)
       lo = lwb(get_box(mf(1),i))
       hi = upb(get_box(mf(1),i))
       select case (dm)
       case (2)
          call nodal_sum_2d(dp(:,:,1,1),ng,sum_proc,vol_proc,lo,hi)
       case (3)
          call nodal_sum_3d(dp(:,:,:,1),ng,sum_proc,vol_proc,lo,hi)
       end select
    end do

    call parallel_reduce(sum, sum_proc, MPI_SUM)
    call parallel_reduce(vol, vol_proc, MPI_SUM)

  end subroutine nodal_sum

  subroutine nodal_sum_2d(mf,ng,sum,vol,lo,hi)

    integer         :: lo(:),hi(:),ng
    real(kind=dp_t) :: sum,vol
    real(kind=dp_t) :: mf(lo(1)-ng:,lo(2)-ng:)

    integer :: i,j

    ! interior
    do j=lo(2)+1,hi(2)
       do i=lo(1)+1,hi(1)
          sum = sum + mf(i,j)
          vol = vol + 1.d0
       end do
    end do

    ! x-faces
    do j=lo(2)+1,hi(2)
       sum = sum + 0.5d0*mf(lo(1)  ,j)
       sum = sum + 0.5d0*mf(hi(1)+1,j)
       vol = vol + 1.d0
    end do

    ! y faces
    do i=lo(1)+1,hi(1)
       sum = sum + 0.5d0*mf(i,lo(2)  )
       sum = sum + 0.5d0*mf(i,hi(2)+1)
       vol = vol + 1.d0
    end do

    ! corners
    sum = sum + 0.25d0*mf(lo(1)  ,lo(2)  )
    sum = sum + 0.25d0*mf(hi(1)+1,lo(2)  )
    sum = sum + 0.25d0*mf(lo(1)  ,hi(2)+1)
    sum = sum + 0.25d0*mf(hi(1)+1,hi(2)+1)
    vol = vol + 1.d0

  end subroutine nodal_sum_2d

  subroutine nodal_sum_3d(mf,ng,sum,vol,lo,hi)

    integer         :: lo(:),hi(:),ng
    real(kind=dp_t) :: sum,vol
    real(kind=dp_t) :: mf(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)

    integer :: i,j,k

    ! interior
    do k=lo(3)+1,hi(3)
       do j=lo(2)+1,hi(2)
          do i=lo(1)+1,hi(1)
             sum = sum + mf(i,j,k)
             vol = vol + 1.d0
          end do
       end do
    end do

    ! x-faces
    do k=lo(3)+1,hi(3)
       do j=lo(2)+1,hi(2)
          sum = sum + 0.5d0*mf(lo(1)  ,j,k)
          sum = sum + 0.5d0*mf(hi(1)+1,j,k)
          vol = vol + 1.d0
       end do
    end do

    ! y-faces
    do k=lo(3)+1,hi(3)
       do i=lo(1)+1,hi(1)
          sum = sum + 0.5d0*mf(i,lo(2)  ,k)
          sum = sum + 0.5d0*mf(i,hi(2)+1,k)
          vol = vol + 1.d0
       end do
    end do

    ! z-faces
    do j=lo(2)+1,hi(2)
       do i=lo(1)+1,hi(1)
          sum = sum + 0.5d0*mf(i,j,lo(3)  )
          sum = sum + 0.5d0*mf(i,j,hi(3)+1)
          vol = vol + 1.d0
       end do
    end do

    ! xy edges
    do k=lo(3)+1,hi(3)
       sum = sum + 0.25d0*mf(lo(1)  ,lo(2)  ,k)
       sum = sum + 0.25d0*mf(lo(1)  ,hi(2)+1,k)
       sum = sum + 0.25d0*mf(hi(1)+1,lo(2)  ,k)
       sum = sum + 0.25d0*mf(hi(1)+1,hi(2)+1,k)
       vol = vol + 1.d0
    end do

    ! xz edges
    do j=lo(2)+1,hi(2)
       sum = sum + 0.25d0*mf(lo(1)  ,j,lo(3)  )
       sum = sum + 0.25d0*mf(lo(1)  ,j,hi(3)+1)
       sum = sum + 0.25d0*mf(hi(1)+1,j,lo(3)  )
       sum = sum + 0.25d0*mf(hi(1)+1,j,hi(3)+1)
       vol = vol + 1.d0
    end do

    ! yz edges
    do i=lo(1)+1,hi(1)
       sum = sum + 0.25d0*mf(i,lo(2)  ,lo(3)  )
       sum = sum + 0.25d0*mf(i,lo(2)  ,hi(3)+1)
       sum = sum + 0.25d0*mf(i,hi(2)+1,lo(3)  )
       sum = sum + 0.25d0*mf(i,hi(2)+1,hi(3)+1)
       vol = vol + 1.d0
    end do

    ! corners
    sum = sum + 0.125d0*mf(lo(1)  ,lo(2)  ,lo(3)  )
    sum = sum + 0.125d0*mf(lo(1)  ,lo(2)  ,hi(3)+1)
    sum = sum + 0.125d0*mf(lo(1)  ,hi(2)+1,lo(3)  )
    sum = sum + 0.125d0*mf(lo(1)  ,hi(2)+1,hi(3)+1)
    sum = sum + 0.125d0*mf(hi(1)+1,lo(2)  ,lo(3)  )
    sum = sum + 0.125d0*mf(hi(1)+1,lo(2)  ,hi(3)+1)
    sum = sum + 0.125d0*mf(hi(1)+1,hi(2)+1,lo(3)  )
    sum = sum + 0.125d0*mf(hi(1)+1,hi(2)+1,hi(3)+1)
    vol = vol + 1.d0

  end subroutine nodal_sum_3d

end module nodal_sum_module
