subroutine compute_flux(phi, ng_p, fluxx, fluxy, ng_f, lo, hi, dx)

  implicit none

  integer lo(2),hi(2),ng_p,ng_f
  double precision   phi(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision fluxx(lo(1)-ng_f:hi(1)+ng_f+1,lo(2)-ng_f:hi(2)+ng_f)
  double precision fluxy(lo(1)-ng_f:hi(1)+ng_f,lo(2)-ng_f:hi(2)+ng_f+1)
  double precision dx

  ! local variables
  integer i,j

  ! x-fluxes
  do j=lo(2),hi(2)
     do i=lo(1),hi(1)+1
        fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / dx
     end do
  end do

  ! y-fluxes
  do j=lo(2),hi(2)+1
     do i=lo(1),hi(1)
        fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / dx
     end do
  end do

end subroutine compute_flux

subroutine update_phi(phiold, phinew, ng_p, fluxx, fluxy, ng_f, lo, hi, dx, dt)

  integer          :: lo(2), hi(2), ng_p, ng_f
  double precision :: phiold(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision :: phinew(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision ::  fluxx(lo(1)-ng_f:hi(1)+ng_f+1,lo(2)-ng_f:hi(2)+ng_f)
  double precision ::  fluxy(lo(1)-ng_f:hi(1)+ng_f,lo(2)-ng_f:hi(2)+ng_f+1)
  double precision :: dx, dt

  ! local variables
  integer i,j

  do j=lo(2),hi(2)
     do i=lo(1),hi(1)

        phinew(i,j) = phiold(i,j) + dt * &
             ( fluxx(i+1,j)-fluxx(i,j) + fluxy(i,j+1)-fluxy(i,j) ) / dx

     end do
  end do

end subroutine update_phi
