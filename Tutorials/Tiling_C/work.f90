
subroutine work(flag, lo, hi, a, a_l1, a_l2, a_l3, a_h1, a_h2, a_h3) bind(C, name="work")

  integer, intent(in) :: flag, lo(3), hi(3), a_l1, a_l2, a_l3, a_h1, a_h2, a_h3
  double precision :: a(a_l1:a_h1, a_l2:a_h2, a_l3:a_h3)
  
  integer :: i, j, k

  if (flag .eq. 0) then
     do k       = lo(3), hi(3)
        do j    = lo(2), hi(2)
           do i = lo(1), hi(1)
              a(i,j,k) = exp(-dble(i*i+j*j+k*k)/512.d0)
           end do
        end do
     end do
  else if (flag .eq. 1) then
     do k       = lo(3), hi(3)
        do j    = lo(2), hi(2)
           do i = lo(1), hi(1)
              a(i,j,k) = log(a(i,j,k))
           end do
        end do
     end do
  else if (flag .eq. 2) then
     do k       = lo(3), hi(3)
        do j    = lo(2), hi(2)
           do i = lo(1), hi(1)
              a(i,j,k) = exp(a(i,j,k))
           end do
        end do
     end do
  end if

end subroutine work
