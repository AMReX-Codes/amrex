subroutine SetIC(mf, lo, hi, dptr) bind(C, name="FSetIC")

  use amrex_error_module
  use rhs_mod
  use ode_params
  use, intrinsic :: iso_c_binding

  implicit none

  integer, intent(in) :: lo(3),hi(3)
  real(c_double), intent(inout) :: mf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
  real(c_double), intent(inout) :: dptr((hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)*(hi(2)-lo(2)+1))
!  real(c_double), intent(inout) :: dptr(0:((hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)*(hi(2)-lo(2)+1)-1))

  ! local variables
  integer i,j,k,n 
  integer(c_long) :: offset, tile_size(3)
  tile_size = hi - lo + 1
  n=1
  print*, tile_size(1)*tile_size(2)*tile_size(3)
  print*, (hi(3)-lo(3)+1)*(hi(2)-lo(2)+1)*(hi(1)-lo(1)+1)
  print*, lo(1), hi(1)
  print*, lo(2), hi(2)
  print*, lo(3), hi(3)

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        offset = tile_size(1) * (j - lo(2)) + tile_size(1) * tile_size(2) * (k - lo(3)) + &
                      tile_size(1) * tile_size(2) * tile_size(3) * (n - 1) + 1 - lo(1)
        do i=lo(1),hi(1)
!           mf(i,j,k) = dptr(offset+1)
           ! Set initial conditions for the ODE for this cell.
           dptr(offset+i) = real(i+j+k, c_double)
        end do
     end do
  end do

end subroutine SetIC

subroutine SetIC_mfab(mf, lo, hi) bind(C, name="FSetIC_mfab")

  use amrex_error_module
  use rhs_mod
  use ode_params
  use, intrinsic :: iso_c_binding

  implicit none

  integer, intent(in) :: lo(3),hi(3)
  real(c_double), intent(inout) :: mf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

  ! local variables
  integer i,j,k,n 
  integer(c_long) :: offset, tile_size(3)
  tile_size = hi - lo + 1
  n=1

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           ! Set initial conditions for the ODE for this cell.
           mf(i,j,k) = real(i+j+k, c_double)
        end do
     end do
  end do

end subroutine SetIC_mfab

subroutine SetSol(mf, lo, hi, dptr) bind(C, name="FSetSol")

  use amrex_error_module
  use rhs_mod
  use ode_params
  use, intrinsic :: iso_c_binding

  implicit none

  integer, intent(in) :: lo(3),hi(3)
  real(c_double), intent(inout) :: mf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
!  real(c_double), intent(inout) :: dptr(0:((hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)*(hi(2)-lo(2)+1)-1))
  real(c_double), intent(inout) :: dptr((hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)*(hi(2)-lo(2)+1))

  ! local variables
  integer i,j,k,n 
  integer(c_long) :: offset, tile_size(3)
  tile_size = hi - lo + 1
  n=1
  print*, tile_size(1)*tile_size(2)*tile_size(3)
  print*, (hi(3)-lo(3)+1)*(hi(2)-lo(2)+1)*(hi(1)-lo(1)+1)
  print*, lo(1), hi(1)
  print*, lo(2), hi(2)
  print*, lo(3), hi(3)

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        offset = tile_size(1) * (j - lo(2)) + tile_size(1) * tile_size(2) * (k - lo(3)) + &
                      tile_size(1) * tile_size(2) * tile_size(3) * (n - 1) + 1 - lo(1)
        do i=lo(1),hi(1)

           ! Set solution for the ODE for this cell.
           mf(i,j,k) = dptr(offset+i)
!           print*, dptr(offset+1)
        end do
     end do
  end do

end subroutine SetSol
