subroutine init_phi(lo, hi, &
     phi, p_l1, p_l2, p_l3, p_h1, p_h2, p_h3, &
     ncomp, dx, prob_lo, prob_hi)

  use my_module
  use mempool_module

  implicit none

  integer, intent(in) :: lo(3), hi(3), ncomp
  integer, intent(in) :: p_l1, p_l2, p_l3, p_h1, p_h2, p_h3
  double precision, intent(inout) :: phi(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3,ncomp)
  double precision :: dx(3), prob_lo(3), prob_hi(3) 

  integer          :: i,j,k,n
  double precision :: x,y,z,r2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  double precision, allocatable :: a(:,:)
  double precision, pointer :: b(:,:)
  double precision, pointer :: c(:,:)
  integer :: mylo(2), myhi(2)

  mylo = 8
  myhi = 15

  allocate(a(mylo(1):myhi(1),mylo(2):myhi(2)))
  allocate(b(mylo(1):myhi(1),mylo(2):myhi(2)))
  call bl_allocate(c, mylo(1),myhi(1),mylo(2),myhi(2))

  a = 1.d0
  b = 1.d0
  c = 1.d0

  print *, '***** For allocatable *****'

  print *, '---- in main ----'

  i = (mylo(1)+myhi(1))/2
  do j=mylo(2),myhi(2)
     print *, loc(a(i,j))
  end do

  call f(mylo, myhi, a)
  call g(mylo, myhi, a)

  print *, '***** For pointer *****'

  print *, '---- in main ----'

  i = (mylo(1)+myhi(1))/2
  do j=mylo(2),myhi(2)
     print *, loc(b(i,j))
  end do

  call f(mylo, myhi, b)
  call g(mylo, myhi, b)

  print *, '***** For bl pointer *****'

  print *, '---- in main ----'

  i = (mylo(1)+myhi(1))/2
  do j=mylo(2),myhi(2)
     print *, loc(c(i,j))
  end do

  call f(mylo, myhi, c)
  call g(mylo, myhi, c)

  deallocate(a,b)
  call bl_deallocate(c)

  stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


