
module amrex_eb2_2d_moudle

  use amrex_error_module
  use amrex_fort_module
  use amrex_constants_module, only : zero, one, half, fourth, sixth, eighth
  use amrex_ebcellflag_module, only : regular_cell => regular, covered_cell => covered, &
       is_regular_cell, is_single_valued_cell, is_covered_cell, get_cell_type, get_neighbor_cells, &
       set_regular_cell, set_single_valued_cell, set_covered_cell, &
       set_neighbor, clear_neighbor, is_neighbor
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
       apx, xlo, xhi, apy, ylo, yhi, cflag, flo, fhi) &
       bind(c,name='amrex_eb2_build_cellflag_from_ap')
    integer, dimension(2), intent(in) :: lo, hi, xlo, xhi, ylo, yhi, flo, fhi
    real(amrex_real), intent(in   ) ::  apx (xlo(1):xhi(1),xlo(2):xhi(2))
    real(amrex_real), intent(in   ) ::  apy (ylo(1):yhi(1),ylo(2):yhi(2))
    integer         , intent(inout) :: cflag(flo(1):fhi(1),flo(2):fhi(2))

    integer :: i,j, flg

    ! By default, all neighbors are already set.
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          flg = cflag(i,j)

          if (apx(i  ,j  ).eq.zero) flg = clear_neighbor(flg, -1,  0)
          if (apx(i+1,j  ).eq.zero) flg = clear_neighbor(flg,  1,  0)
          if (apy(i  ,j  ).eq.zero) flg = clear_neighbor(flg,  0, -1)
          if (apy(i  ,j+1).eq.zero) flg = clear_neighbor(flg,  0,  1)

          if (apx(i,j).ne.zero .and. apy(i-1,j).ne.zero) then
          else if (apx(i,j-1).ne.zero .and. apy(i,j).ne.zero) then
          else
             flg = clear_neighbor(flg,-1,-1)
          end if

          if (apx(i+1,j).ne.zero .and. apy(i+1,j).ne.zero) then
          else if (apx(i+1,j-1).ne.zero .and. apy(i,j).ne.zero) then
          else
             flg = clear_neighbor(flg,1,-1)
          end if

          if (apx(i,j).ne.zero .and. apy(i-1,j+1).ne.zero) then
          else if (apx(i,j+1).ne.zero .and. apy(i,j+1).ne.zero) then
          else
             flg = clear_neighbor(flg,-1,1)
          end if

          if (apx(i+1,j).ne.zero .and. apy(i+1,j+1).ne.zero) then
          else if (apx(i+1,j+1).ne.zero .and. apy(i,j+1).ne.zero) then
          else
             flg = clear_neighbor(flg,1,1)
          end if

          cflag(i,j) = flg
       end do
    end do
  end subroutine amrex_eb2_build_cellflag_from_ap

end module amrex_eb2_2d_moudle
