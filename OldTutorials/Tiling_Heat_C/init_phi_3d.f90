subroutine init_phi(lo, hi, &
     phi, p_l1, p_l2, p_l3, p_h1, p_h2, p_h3, &
     ncomp, dx, prob_lo, prob_hi) bind(C, name="init_phi")

  implicit none

  integer, intent(in) :: lo(3), hi(3), ncomp
  integer, intent(in) :: p_l1, p_l2, p_l3, p_h1, p_h2, p_h3
  double precision, intent(inout) :: phi(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3,ncomp)
  double precision :: dx(3), prob_lo(3), prob_hi(3) 

  integer          :: i,j,k,n
  double precision :: x,y,z,r2

  do n = 1, ncomp
     do k = lo(3), hi(3)
        z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
        do j = lo(2), hi(2)
           y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
           do i = lo(1), hi(1)
              x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
              
              r2 = ((x-0.25d0)**2 + (y-0.25d0)**2 + (z-0.25d0)**2) * 100.d0
              phi(i,j,k,n) = 1.d0 + exp(-r2)
           end do
        end do
     end do
  end do

end subroutine init_phi
