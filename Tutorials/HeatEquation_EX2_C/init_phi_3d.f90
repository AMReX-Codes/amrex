subroutine init_phi(phi, lo, hi, ng, dx, prob_lo, prob_hi)

  implicit none

  integer          :: lo(3), hi(3), ng
  double precision :: phi(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
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

           r2 = ((x-0.25d0)**2 + (y-0.25d0)**2 + (z-0.25d0)**2) / 0.01d0
           phi(i,j,k) = 1.d0 + exp(-r2)

        end do
     end do
  end do

end subroutine init_phi
