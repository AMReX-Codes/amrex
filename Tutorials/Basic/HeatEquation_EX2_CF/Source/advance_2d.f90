
subroutine compute_flux (lo, hi, domlo, domhi, phi, philo, phihi, &
                         fluxx, fxlo, fxhi, fluxy, fylo, fyhi, &
                         dx, bc) bind(C, name="compute_flux")

  use amrex_fort_module, only : amrex_real
  use amrex_bc_types_module
  implicit none

  integer lo(2), hi(2), domlo(2), domhi(2)
  integer philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  integer, intent(in)             :: bc(2,2,1) ! (dim,lohi,ncomp)

  ! local variables
  integer i,j

  ! x-fluxes
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)+1
     fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / dx(1)
  end do
  end do

  ! y-fluxes
  do j = lo(2), hi(2)+1
  do i = lo(1), hi(1)
     fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / dx(2)
  end do
  end do

  ! lo-x boundary, ghost cell contains value on boundary
  if (domlo(1) .eq. lo(1) .and. &
       (bc(1,1,1) .eq. amrex_bc_foextrap .or. bc(1,1,1) .eq. amrex_bc_ext_dir) ) then
     i = lo(1)
     do j = lo(2), hi(2)
        fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / (0.5d0*dx(1))
     end do
  end if

  ! hi-x boundary, ghost cell contains value on boundary
  if (domhi(1) .eq. hi(1) .and. &
       (bc(1,2,1) .eq. amrex_bc_foextrap .or. bc(1,2,1) .eq. amrex_bc_ext_dir) ) then
     i = hi(1)+1
     do j = lo(2), hi(2)
        fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / (0.5d0*dx(1))
     end do
  end if

  ! lo-y boundary, ghost cell contains value on boundary
  if (domlo(2) .eq. lo(2) .and. &
       (bc(2,1,1) .eq. amrex_bc_foextrap .or. bc(2,1,1) .eq. amrex_bc_ext_dir) ) then
     j = lo(2)
     do i = lo(1), hi(1)
        fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / (0.5d0*dx(2))
     end do
  end if

  ! lo-y boundary, ghost cell contains value on boundary
  if (domhi(2) .eq. hi(2) .and. &
       (bc(2,2,1) .eq. amrex_bc_foextrap .or. bc(2,2,1) .eq. amrex_bc_ext_dir) ) then
     j = hi(2)+1
     do i = lo(1), hi(1)
        fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / (0.5d0*dx(2))
     end do
  end if

end subroutine compute_flux


subroutine update_phi (lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
                       fluxx, fxlo, fxhi, fluxy, fylo, fyhi, &
                       dx, dt) bind(C, name="update_phi")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), polo(2), pohi(2), pnlo(2), pnhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
  real(amrex_real), intent(in)    :: phiold(polo(1):pohi(1),polo(2):pohi(2))
  real(amrex_real), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2))
  real(amrex_real), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
  real(amrex_real), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  real(amrex_real), intent(in)    :: dt

  ! local variables
  integer i,j
  real(amrex_real) :: dtdx(2)

  dtdx = dt/dx

  do j = lo(2), hi(2)
  do i = lo(1), hi(1)

     phinew(i,j) = phiold(i,j) &
          + dtdx(1) * (fluxx(i+1,j  ) - fluxx(i,j)) &
          + dtdx(2) * (fluxy(i  ,j+1) - fluxy(i,j))

  end do
  end do

end subroutine update_phi
