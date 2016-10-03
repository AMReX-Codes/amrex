
subroutine advance(U, Unew, lo, hi, Ncomp, ng, dx, dt) bind(C, name="advance")

  implicit none

  integer lo(3),hi(3),Ncomp,ng
  double precision    U(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng,Ncomp)
  double precision Unew(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng,Ncomp)
  double precision dx,dt

  ! Local variables
  integer ndo
  double precision, allocatable :: dU(:,:,:,:)
  double precision a

  double precision, parameter :: third     =  1.d0/3.d0
  double precision, parameter :: twothirds =  2.d0/3.d0
  double precision, parameter :: fac1      = -1.d0/12.d0
  double precision, parameter :: fac2      =  4.d0/3.d0
  double precision, parameter :: fac3      = -7.5d0

  integer          :: i,j,k

  allocate(dU(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,Ncomp))

  if (ng < 6) then
     print *,'NOT ENOUGH GHOST CELLS IN ADVANCE ',ng
     stop
  end if

  ! First call to get RHS = Lap(U)
  ndo = 4

  do k = lo(3)-ndo,hi(3)+ndo
     do j = lo(2)-ndo,hi(2)+ndo
        do i = lo(1)-ndo,hi(1)+ndo

           dU(i,j,k,1) = U(i,j,k,2)

           dU(i,j,k,2) = (  fac1 * ( U(i-2,j,k,1) + U(i+2,j,k,1) + &
                                     U(i,j-2,k,1) + U(i,j+2,k,1) + &
                                     U(i,j,k-2,1) + U(i,j,k+2,1) ) &
                          + fac2 * ( U(i-1,j,k,1) + U(i+1,j,k,1) + &
                                     U(i,j-1,k,1) + U(i,j+1,k,1) + &
                                     U(i,j,k-1,1) + U(i,j,k+1,1) ) &
                          + fac3 *   U(i,j,k,1)  ) / (dx*dx)
        end do
     end do
  end do

  ! First update
  ! This Unew lives at t^{n+1}
  Unew(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2) = & 
               U(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2) &
       + dt * dU(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2)

  ! Second call to get RHS = Lap(Unew)
  ndo = 2

  do k = lo(3)-ndo,hi(3)+ndo
     do j = lo(2)-ndo,hi(2)+ndo
        do i = lo(1)-ndo,hi(1)+ndo

           dU(i,j,k,1) = Unew(i,j,k,2)

           dU(i,j,k,2) = (  fac1 * ( Unew(i-2,j,k,1) + Unew(i+2,j,k,1) + &
                                     Unew(i,j-2,k,1) + Unew(i,j+2,k,1) + &
                                     Unew(i,j,k-2,1) + Unew(i,j,k+2,1) ) &
                          + fac2 * ( Unew(i-1,j,k,1) + Unew(i+1,j,k,1) + &
                                     Unew(i,j-1,k,1) + Unew(i,j+1,k,1) + &
                                     Unew(i,j,k-1,1) + Unew(i,j,k+1,1) ) &
                          + fac3 *   Unew(i,j,k,1)  ) / (dx*dx)
        end do
     end do
  end do

  ! Second update
  ! This Unew lives at t^{n+1/2}
  Unew(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2) = & 
              .75d0 *    U(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2) &
            + .25d0 * Unew(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2) &
       + dt * .25d0 *   dU(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2) 

  ! Third call to get RHS = Lap(Unew)
  ndo = 0

  do k = lo(3)-ndo,hi(3)+ndo
     do j = lo(2)-ndo,hi(2)+ndo
        do i = lo(1)-ndo,hi(1)+ndo

           dU(i,j,k,1) = Unew(i,j,k,2)

           dU(i,j,k,2) = (  fac1 * ( Unew(i-2,j,k,1) + Unew(i+2,j,k,1) + &
                                     Unew(i,j-2,k,1) + Unew(i,j+2,k,1) + &
                                     Unew(i,j,k-2,1) + Unew(i,j,k+2,1) ) &
                          + fac2 * ( Unew(i-1,j,k,1) + Unew(i+1,j,k,1) + &
                                     Unew(i,j-1,k,1) + Unew(i,j+1,k,1) + &
                                     Unew(i,j,k-1,1) + Unew(i,j,k+1,1) ) &
                          + fac3 *   Unew(i,j,k,1)  ) / (dx*dx)

        end do
     end do
  end do

  ! Third update
  ! This Unew lives at t^{n+1}
  Unew(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2) =        &
              third     *    U(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2)   &
            + twothirds * Unew(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2)   &
       + dt * twothirds *   dU(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2)

  deallocate(dU)

end subroutine advance
