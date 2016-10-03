subroutine init_data(U, lo, hi, Ncomp, ng,  &
                     dx, prob_lo, prob_hi) bind(C, name="init_data")

  implicit none

  integer          :: lo(2), hi(2), Ncomp, ng
  double precision :: U(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, Ncomp)
  double precision :: dx(2) 
  double precision :: prob_lo(2) 
  double precision :: prob_hi(2) 

  integer          :: i,j
  double precision :: x,y,r2

  do j = lo(2), hi(2)
     y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
     do i = lo(1), hi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
        r2 = (x*x + y*y) / 0.01

        U(i,j,1) = 0.d0
        U(i,j,2) = exp(-r2)

     end do
  end do

end subroutine init_data
