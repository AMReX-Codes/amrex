subroutine shift_E(E, lo, hi, ng, N) bind(C, name="shift_E")

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

end subroutine shift_E
