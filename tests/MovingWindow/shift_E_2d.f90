subroutine shift_x(E, lo, hi, ng, N) bind(C, name="shift_x")

  implicit none

  integer          :: lo(2), hi(2), ng, N
  double precision :: E(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)

  integer          :: i,j

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        do n = 0, N - 1
           E(i + n, j) = E(i + n + 1, j)
        end do
     end do
  end do

end subroutine shift_x


subroutine shift_y(E, lo, hi, ng, N) bind(C, name="shift_y")

  implicit none

  integer          :: lo(2), hi(2), ng, N
  double precision :: E(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)

  integer          :: i,j

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        do n = 0, N - 1
           E(i, j + n) = E(i, j + n + 1)
        end do
     end do
  end do

end subroutine shift_y
