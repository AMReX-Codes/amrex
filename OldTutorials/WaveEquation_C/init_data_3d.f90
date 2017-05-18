
subroutine init_data(U, lo, hi, Ncomp, ng,  &
                     dx, prob_lo, prob_hi) bind(C, name="init_data")

  implicit none

  integer          :: lo(3), hi(3), Ncomp, ng
  double precision :: U(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, Ncomp)
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
           r2 = (x*x + y*y + z*z) / 0.01

           U(i,j,k,1) = 0.d0
           U(i,j,k,2) = exp(-r2)

        end do
     end do
  end do

end subroutine init_data
