subroutine compute_flux(phi, ng_p, fluxx, fluxy, fluxz, ng_f, lo, hi, dx)

  implicit none

  integer lo(3),hi(3),ng_p,ng_f
  double precision   phi(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision fluxx(lo(1)-ng_f:hi(1)+ng_f+1,lo(2)-ng_f:hi(2)+ng_f,lo(3)-ng_f:hi(3)+ng_f)
  double precision fluxy(lo(1)-ng_f:hi(1)+ng_f,lo(2)-ng_f:hi(2)+ng_f+1,lo(3)-ng_f:hi(3)+ng_f)
  double precision fluxz(lo(1)-ng_f:hi(1)+ng_f,lo(2)-ng_f:hi(2)+ng_f,lo(3)-ng_f:hi(3)+ng_f+1)
  double precision dx
  
  ! local variables
  integer i,j,k

  ! x-fluxes
  !$omp parallel do private(i,j,k)
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)+1
           fluxx(i,j,k) = ( phi(i,j,k) - phi(i-1,j,k) ) / dx
        end do
     end do
  end do
  !$omp end parallel do

  ! y-fluxes
  !$omp parallel do private(i,j,k)
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)+1
        do i=lo(1),hi(1)
           fluxy(i,j,k) = ( phi(i,j,k) - phi(i,j-1,k) ) / dx
        end do
     end do
  end do
  !$omp end parallel do

  ! z-fluxes
  !$omp parallel do private(i,j,k)
  do k=lo(3),hi(3)+1
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           fluxz(i,j,k) = ( phi(i,j,k) - phi(i,j,k-1) ) / dx
        end do
     end do
  end do
  !$omp end parallel do

end subroutine compute_flux

subroutine update_phi(phiold, phinew, ng_p, fluxx, fluxy, fluxz, ng_f, lo, hi, dx, dt)

  integer          :: lo(3), hi(3), ng_p, ng_f
  double precision :: phiold(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision :: phinew(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision ::  fluxx(lo(1)-ng_f:hi(1)+ng_f+1,lo(2)-ng_f:hi(2)+ng_f,lo(3)-ng_f:hi(3)+ng_f)
  double precision ::  fluxy(lo(1)-ng_f:hi(1)+ng_f,lo(2)-ng_f:hi(2)+ng_f+1,lo(3)-ng_f:hi(3)+ng_f)
  double precision ::  fluxz(lo(1)-ng_f:hi(1)+ng_f,lo(2)-ng_f:hi(2)+ng_f,lo(3)-ng_f:hi(3)+ng_f+1)
  double precision :: dx, dt

  ! local variables
  integer i,j,k

  !$omp parallel do private(i,j,k)
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           phinew(i,j,k) = phiold(i,j,k) + dt * &
                ( fluxx(i+1,j,k)-fluxx(i,j,k) &
                +fluxy(i,j+1,k)-fluxy(i,j,k) &
                +fluxz(i,j,k+1)-fluxz(i,j,k) ) / dx

        end do
     end do
  end do
  !$omp end parallel do

end subroutine update_phi
