subroutine work_on_data(data, ng, nc, lo, hi) bind(C, name="work_on_data")

  implicit none

  integer          :: lo(3), hi(3), ng, nc
  double precision :: data(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,1:nc)

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

end subroutine work_on_data
