subroutine advance_phi(lo, hi, &
     pold, po_l1, po_l2, po_l3, po_h1, po_h2, po_h3, &
     pnew, pn_l1, pn_l2, pn_l3, pn_h1, pn_h2, pn_h3, &
     ncomp, dx, dt) bind(C, name="advance_phi")

  implicit none

  integer, intent(in) :: lo(3),hi(3), ncomp
  integer, intent(in) :: po_l1, po_l2, po_l3, po_h1, po_h2, po_h3
  integer, intent(in) :: pn_l1, pn_l2, pn_l3, pn_h1, pn_h2, pn_h3
  double precision, intent(in   ) :: pold(po_l1:po_h1,po_l2:po_h2,po_l3:po_h3,ncomp)
  double precision, intent(inout) :: pnew(pn_l1:pn_h1,pn_l2:pn_h2,pn_l3:pn_h3,ncomp)
  double precision, intent(in   ) :: dx, dt
  
  integer :: i,j,k,n
  double precision :: dxinv, dtdx
  double precision :: fx(lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3))
  double precision :: fy(lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3))
  double precision :: fz(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1)

  dxinv = 1.d0/dx
  dtdx = dt*dxinv

  do n = 1, ncomp
     ! x-fluxes
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)+1
              fx(i,j,k) = ( pold(i-1,j,k,n) - pold(i,j,k,n) ) * dxinv
           end do
        end do
     end do

     ! y-fluxes
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)+1
           do i=lo(1),hi(1)
              fy(i,j,k) = ( pold(i,j-1,k,n) - pold(i,j,k,n) ) * dxinv
           end do
        end do
     end do
     
     ! z-fluxes
     do k=lo(3),hi(3)+1
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)
              fz(i,j,k) = ( pold(i,j,k-1,n) - pold(i,j,k,n) ) * dxinv
           end do
        end do
     end do
     
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)
              
              pnew(i,j,k,n) = pold(i,j,k,n) - dtdx * &
                   ( fx(i+1,j,k)-fx(i,j,k) &
                   + fy(i,j+1,k)-fy(i,j,k) &
                   + fz(i,j,k+1)-fz(i,j,k) )

           end do
        end do
     end do

  end do

end subroutine advance_phi


subroutine advance_phi2(lo, hi, &
     pold, po_l1, po_l2, po_l3, po_h1, po_h2, po_h3, &
     pnew, pn_l1, pn_l2, pn_l3, pn_h1, pn_h2, pn_h3, &
     ncomp, dx, dt) bind(C, name="advance_phi2")

  implicit none

  integer, intent(in) :: lo(3),hi(3), ncomp
  integer, intent(in) :: po_l1, po_l2, po_l3, po_h1, po_h2, po_h3
  integer, intent(in) :: pn_l1, pn_l2, pn_l3, pn_h1, pn_h2, pn_h3
  double precision, intent(in   ) :: pold(po_l1:po_h1,po_l2:po_h2,po_l3:po_h3,ncomp)
  double precision, intent(inout) :: pnew(pn_l1:pn_h1,pn_l2:pn_h2,pn_l3:pn_h3,ncomp)
  double precision, intent(in   ) :: dx, dt
  
  integer :: i,j,k,n
  double precision :: dxinv, dtdx
  double precision :: fx(lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3))
  double precision :: fy(lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3))
  double precision :: fz(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1)

  dxinv = 1.d0/dx
  dtdx = dt*dxinv

  !$omp parallel private(i,j,k,n)

  do n = 1, ncomp
     ! x-fluxes
     !$omp do collapse(2)
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)+1
              fx(i,j,k) = ( pold(i-1,j,k,n) - pold(i,j,k,n) ) * dxinv
           end do
        end do
     end do
     !$omp end do nowait

     ! y-fluxes
     !$omp do collapse(2)
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)+1
           do i=lo(1),hi(1)
              fy(i,j,k) = ( pold(i,j-1,k,n) - pold(i,j,k,n) ) * dxinv
           end do
        end do
     end do
     !$omp end do nowait
     
     ! z-fluxes
     !$omp do collapse(2)
     do k=lo(3),hi(3)+1
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)
              fz(i,j,k) = ( pold(i,j,k-1,n) - pold(i,j,k,n) ) * dxinv
           end do
        end do
     end do
     !$omp end do
     
     !$omp do collapse(2)
     do k=lo(3),hi(3)
        do j=lo(2),hi(2)
           do i=lo(1),hi(1)
              
              pnew(i,j,k,n) = pold(i,j,k,n) - dtdx * &
                   ( fx(i+1,j,k)-fx(i,j,k) &
                   + fy(i,j+1,k)-fy(i,j,k) &
                   + fz(i,j,k+1)-fz(i,j,k) )

           end do
        end do
     end do
     !$omp end do

  end do

  !$omp end parallel

end subroutine advance_phi2
