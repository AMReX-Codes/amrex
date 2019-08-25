
module amrex_eb2_3d_module

  use amrex_error_module
  use amrex_fort_module
  use amrex_constants_module, only : zero, one, two, three, four, five, six, seven, eight,&
       nine, ten, eleven, fifteen, sixteen, half, third, fourth, sixth, eighth, twelfth
  use amrex_ebcellflag_module, only : regular_cell => regular, covered_cell => covered, &
       is_regular_cell, is_single_valued_cell, is_covered_cell, get_cell_type, get_neighbor_cells, &
       set_regular_cell, set_single_valued_cell, set_covered_cell, &
       set_neighbor, clear_neighbor, clear_allneighbors, is_neighbor
  implicit none

  integer, parameter :: regular = 0
  integer, parameter :: covered = 1
  integer, parameter :: irregular = 2
  integer, parameter :: unknown = 3

  real(amrex_real), private, parameter :: small = 1.d-14

  private
  public :: amrex_eb2_build_cellflag_from_ap

contains

  subroutine amrex_eb2_build_cellflag_from_ap (lo, hi, &
       apx, xlo, xhi, apy, ylo, yhi, apz, zlo, zhi, cflag, flo, fhi) &
       bind(c,name='amrex_eb2_build_cellflag_from_ap')
    integer, dimension(3) :: lo, hi, xlo, xhi, ylo, yhi, zlo, zhi, flo, fhi
    real(amrex_real), intent(in   ) ::  apx (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    real(amrex_real), intent(in   ) ::  apy (ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3))
    real(amrex_real), intent(in   ) ::  apz (zlo(1):zhi(1),zlo(2):zhi(2),zlo(3):zhi(3))
    integer         , intent(inout) :: cflag(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

    integer :: i,j,k, flg

    ! By default, all neighbors are already set.
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             cflag(i,j,k) = clear_allneighbors(cflag(i,j,k))
                
             if (.not.is_covered_cell(cflag(i,j,k))) then
                flg = clear_allneighbors(cflag(i,j,k))
                
                if (apx(i  ,j,k).ne.zero) flg = set_neighbor(flg, -1,  0,  0)
                if (apx(i+1,j,k).ne.zero) flg = set_neighbor(flg,  1,  0,  0)
                if (apy(i,j  ,k).ne.zero) flg = set_neighbor(flg,  0, -1,  0)
                if (apy(i,j+1,k).ne.zero) flg = set_neighbor(flg,  0,  1,  0)
                if (apz(i,j,k  ).ne.zero) flg = set_neighbor(flg,  0,  0, -1)
                if (apz(i,j,k+1).ne.zero) flg = set_neighbor(flg,  0,  0,  1)
                
                if ( (apx(i,j,k).ne.zero .and. apy(i-1,j,k).ne.zero) .or. &
                     (apy(i,j,k).ne.zero .and. apx(i,j-1,k).ne.zero) ) then
                   flg = set_neighbor(flg, -1, -1, 0)
                   if (apz(i-1,j-1,k  ).ne.zero) flg = set_neighbor(flg,-1,-1,-1)
                   if (apz(i-1,j-1,k+1).ne.zero) flg = set_neighbor(flg,-1,-1, 1)
                end if
                
                if ( (apx(i+1,j,k).ne.zero .and. apy(i+1,j,k).ne.zero) .or. &
                     (apy(i,j,k).ne.zero .and. apx(i+1,j-1,k).ne.zero) ) then
                   flg = set_neighbor(flg, 1, -1, 0)
                   if (apz(i+1,j-1,k  ).ne.zero) flg = set_neighbor(flg, 1,-1,-1)
                   if (apz(i+1,j-1,k+1).ne.zero) flg = set_neighbor(flg, 1,-1, 1)
                end if
                
                if ( (apx(i,j,k).ne.zero .and. apy(i-1,j+1,k).ne.zero) .or. &
                     (apy(i,j+1,k).ne.zero .and. apx(i,j+1,k).ne.zero) ) then
                   flg = set_neighbor(flg, -1, 1, 0)
                   if (apz(i-1,j+1,k  ).ne.zero) flg = set_neighbor(flg,-1, 1,-1)
                   if (apz(i-1,j+1,k+1).ne.zero) flg = set_neighbor(flg,-1, 1, 1)
                end if
                
                if ( (apx(i+1,j,k).ne.zero .and. apy(i+1,j+1,k).ne.zero) .or. &
                     (apy(i,j+1,k).ne.zero .and. apx(i+1,j+1,k).ne.zero) ) then
                   flg = set_neighbor(flg, 1, 1, 0)
                   if (apz(i+1,j+1,k  ).ne.zero) flg = set_neighbor(flg, 1, 1,-1)
                   if (apz(i+1,j+1,k+1).ne.zero) flg = set_neighbor(flg, 1, 1, 1)
                end if
                
                if ( (apx(i,j,k).ne.zero .and. apz(i-1,j,k).ne.zero) .or. &
                     (apz(i,j,k).ne.zero .and. apx(i,j,k-1).ne.zero) ) then
                   flg = set_neighbor(flg, -1, 0, -1)
                   if (apy(i-1,j  ,k-1).ne.zero) flg = set_neighbor(flg,-1,-1,-1)
                   if (apy(i-1,j+1,k-1).ne.zero) flg = set_neighbor(flg,-1, 1,-1)
                end if
                
                if ( (apx(i+1,j,k).ne.zero .and. apz(i+1,j,k).ne.zero) .or. &
                     (apz(i,j,k).ne.zero .and. apx(i+1,j,k-1).ne.zero) ) then
                   flg = set_neighbor(flg, 1, 0, -1)
                   if (apy(i+1,j  ,k-1).ne.zero) flg = set_neighbor(flg, 1,-1,-1)
                   if (apy(i+1,j+1,k-1).ne.zero) flg = set_neighbor(flg, 1, 1,-1)
                end if
                
                if ( (apx(i,j,k).ne.zero .and. apz(i-1,j,k+1).ne.zero) .or. &
                     (apz(i,j,k+1).ne.zero .and. apx(i,j,k+1).ne.zero) ) then
                   flg = set_neighbor(flg, -1, 0, 1)
                   if (apy(i-1,j  ,k+1).ne.zero) flg = set_neighbor(flg,-1,-1, 1)
                   if (apy(i-1,j+1,k+1).ne.zero) flg = set_neighbor(flg,-1, 1, 1)
                end if
                
                if ( (apx(i+1,j,k).ne.zero .and. apz(i+1,j,k+1).ne.zero) .or. &
                     (apz(i,j,k+1).ne.zero .and. apx(i+1,j,k+1).ne.zero) ) then
                   flg = set_neighbor(flg, 1, 0, 1)
                   if (apy(i+1,j  ,k+1).ne.zero) flg = set_neighbor(flg, 1,-1, 1)
                   if (apy(i+1,j+1,k+1).ne.zero) flg = set_neighbor(flg, 1, 1, 1)
                end if
                
                if ( (apy(i,j,k).ne.zero .and. apz(i,j-1,k).ne.zero) .or. &
                     (apz(i,j,k).ne.zero .and. apy(i,j,k-1).ne.zero) ) then
                   flg = set_neighbor(flg, 0, -1, -1)
                   if (apx(i  ,j-1,k-1).ne.zero) flg = set_neighbor(flg,-1,-1,-1)
                   if (apx(i+1,j-1,k-1).ne.zero) flg = set_neighbor(flg, 1,-1,-1)
                end if
                
                if ( (apy(i,j+1,k).ne.zero .and. apz(i,j+1,k).ne.zero) .or. &
                     (apz(i,j,k).ne.zero .and. apy(i,j+1,k-1).ne.zero) ) then
                   flg = set_neighbor(flg, 0, 1, -1)
                   if (apx(i  ,j+1,k-1).ne.zero) flg = set_neighbor(flg,-1, 1,-1)
                   if (apx(i+1,j+1,k-1).ne.zero) flg = set_neighbor(flg, 1, 1,-1)
                end if
                
                if ( (apy(i,j,k).ne.zero .and. apz(i,j-1,k+1).ne.zero) .or. &
                     (apz(i,j,k+1).ne.zero .and. apy(i,j,k+1).ne.zero) ) then
                   flg = set_neighbor(flg, 0, -1, 1)
                   if (apx(i  ,j-1,k+1).ne.zero) flg = set_neighbor(flg,-1,-1, 1)
                   if (apx(i+1,j-1,k+1).ne.zero) flg = set_neighbor(flg, 1,-1, 1)
                end if
                
                if ( (apy(i,j+1,k).ne.zero .and. apz(i,j+1,k+1).ne.zero) .or. &
                     (apz(i,j,k+1).ne.zero .and. apy(i,j+1,k+1).ne.zero) ) then
                   flg = set_neighbor(flg, 0, 1, 1)
                   if (apx(i  ,j+1,k+1).ne.zero) flg = set_neighbor(flg,-1, 1, 1)
                   if (apx(i+1,j+1,k+1).ne.zero) flg = set_neighbor(flg, 1, 1, 1)
                end if
                   
                cflag(i,j,k) = flg
             end if
          end do
       end do
    end do
  end subroutine amrex_eb2_build_cellflag_from_ap
  
end module amrex_eb2_3d_module
