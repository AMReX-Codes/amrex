subroutine compute_flux(phi, ng_p, fluxx, fluxy, ng_f, lo, hi, &
                        domlo, domhi, bc, dx) bind(C, name="compute_flux")


  implicit none

  ! includes definitions found in BC_TYPES.H
  include 'bc_types.fi'

  integer lo(2),hi(2),domlo(2),domhi(2),bc(2,2),ng_p,ng_f
  double precision   phi(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision fluxx(lo(1)-ng_f:hi(1)+ng_f+1,lo(2)-ng_f:hi(2)+ng_f)
  double precision fluxy(lo(1)-ng_f:hi(1)+ng_f,lo(2)-ng_f:hi(2)+ng_f+1)
  double precision dx

  ! local variables
  integer i,j
  integer nlft, nrgt, nbot, ntop

  ! x-fluxes
  do j=lo(2),hi(2)
     do i=lo(1),hi(1)+1
        fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / dx
     end do
  end do

  ! lo-x boundary conditions
  if (lo(1) .eq. domlo(1)) then
     if (bc(1,1) .eq. EXT_DIR) then
        i = lo(1)
        do j=lo(2),hi(2)
           ! divide by 0.5*dx since the ghost cell value represents
           ! the value at the wall, not the ghost cell-center
           fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / (0.5d0*dx)
        end do
     else if (bc(1,1) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxx(lo(1),lo(2):hi(2)) = 0.d0
    end if
  end if

  ! hi-x boundary conditions
  if (hi(1) .eq. domhi(1)) then
     if (bc(1,2) .eq. EXT_DIR) then
        i = hi(1)+1
        do j=lo(2),hi(2)
           ! divide by 0.5*dx since the ghost cell value represents
           ! the value at the wall, not the ghost cell-center
           fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / (0.5d0*dx)
        end do
     else if (bc(1,2) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxx(hi(1)+1,lo(2):hi(2)) = 0.d0
    end if
  end if

  ! y-fluxes
  do j=lo(2),hi(2)+1
     do i=lo(1),hi(1)
        fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / dx
     end do
  end do

  ! lo-y boundary conditions
  if (lo(2) .eq. domlo(2)) then
     if (bc(2,1) .eq. EXT_DIR) then
        j = lo(2)
        do i=lo(1),hi(1)
           ! divide by 0.5*dx since the ghost cell value represents
           ! the value at the wall, not the ghost cell-center
           fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / (0.5d0*dx)
        end do
     else if (bc(2,1) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxy(lo(1):hi(1),lo(2)) = 0.d0
    end if
  end if

  ! hi-y boundary conditions
  if (hi(2) .eq. domhi(2)) then
     if (bc(2,2) .eq. EXT_DIR) then
        j = hi(2)+1
        do i=lo(1),hi(1)
           ! divide by 0.5*dx since the ghost cell value represents
           ! the value at the wall, not the ghost cell-center
           fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / (0.5d0*dx)
        end do
     else if (bc(2,2) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxy(lo(1):hi(1),hi(2)+1) = 0.d0
    end if
  end if

end subroutine compute_flux

subroutine update_phi(phiold, phinew, ng_p, fluxx, fluxy, ng_f, lo, hi, dx, dt) bind(C, name="update_phi")

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
             ( fluxx(i+1,j)-fluxx(i,j) &
              +fluxy(i,j+1)-fluxy(i,j) ) / dx

     end do
  end do

end subroutine update_phi
