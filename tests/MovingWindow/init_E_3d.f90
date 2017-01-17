subroutine init_E(E, lo, hi, ng, dx, prob_lo, prob_hi) bind(C, name="init_E")

  implicit none

  integer          :: lo(3), hi(3), ng
  double precision :: E(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: dx(3) 
  double precision :: prob_lo(3) 
  double precision :: prob_hi(3) 

  integer          :: i,j,k
  double precision :: x,y,z,r2

  do k = lo(3), hi(3)
     z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
     do j = lo(2), hi(2)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
        do i = lo(1), hi(1)
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

           E(i,j,k) = x

        end do
     end do
  end do

end subroutine init_E
