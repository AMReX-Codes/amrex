subroutine shift_x(E, lo, hi, ng, N) bind(C, name="shift_x")

  implicit none

  integer          :: lo(3), hi(3), ng, N
  double precision :: E(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)

  integer          :: i,j,k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           do n = 0, N - 1
              E(i + n, j, k) = E(i + n + 1, j, k)
           end do
        end do
     end do
  end do

end subroutine shift_x


subroutine shift_y(E, lo, hi, ng, N) bind(C, name="shift_y")

  implicit none

  integer          :: lo(3), hi(3), ng, N
  double precision :: E(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)

  integer          :: i,j,k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           do n = 0, N - 1
              E(i, j + n, k) = E(i, j + n + 1, k)
           end do
        end do
     end do
  end do

end subroutine shift_y


subroutine shift_z(E, lo, hi, ng, N) bind(C, name="shift_z")

  implicit none

  integer          :: lo(3), hi(3), ng, N
  double precision :: E(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)

  integer          :: i,j,k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           do n = 0, N - 1
              E(i, j, k + n) = E(i, j, k + n + 1)
           end do
        end do
     end do
  end do

end subroutine shift_z
