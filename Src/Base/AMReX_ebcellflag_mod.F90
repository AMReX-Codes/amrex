
module amrex_ebcellflag_module
  implicit none
  private
  public :: is_regular_cell, is_single_valued_cell, is_multi_valued_cell, &
       is_covered_cell, get_neighbor_cells, num_neighbor_cells
  
  integer, parameter :: w_type      = 2
  integer, parameter :: w_numvofs   = 3
  integer, parameter :: pos_numvofs = 2;
  integer, parameter :: pos_ngbr(-1:1,-1:1,-1:1) = &
       reshape((/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26/),&
               [3,3,3]) + (w_type+w_numvofs)

  integer, parameter :: regular       = 0
  integer, parameter :: single_valued = 1
  integer, parameter :: multi_valued  = 2
  integer, parameter :: covered       = 3

contains
  
  elemental logical function is_regular_cell (flag)
    integer, intent(in) :: flag
    is_regular_cell = ibits(flag,0,w_type) .eq. regular
  end function is_regular_cell

  elemental logical function is_single_valued_cell (flag)
    integer, intent(in) :: flag
    is_single_valued_cell = ibits(flag,0,w_type) .eq. single_valued
  end function is_single_valued_cell

  elemental logical function is_multi_valued_cell (flag)
    integer, intent(in) :: flag
    is_multi_valued_cell = ibits(flag,0,w_type) .eq. multi_valued
  end function is_multi_valued_cell

  elemental logical function is_covered_cell (flag)
    integer, intent(in) :: flag
    is_covered_cell = ibits(flag,0,w_type) .eq. covered
  end function is_covered_cell

#if (AMREX_SPACEDIM == 2)

  pure subroutine get_neighbor_cells (flag, ngbr)
    integer, intent(in) :: flag
    integer, intent(inout) :: ngbr(-1:1,-1:1)
    integer :: i, j
    do j = -1,1
       do i = -1,1
          if (btest(flag,pos_ngbr(i,j,0))) then
             ngbr(i,j) = 1
          else
             ngbr(i,j) = 0
          end if
       end do
    end do
  end subroutine get_neighbor_cells

  elemental integer function num_neighbor_cells (flag)
    integer, intent(in) :: flag
    integer ngbr(-1:1,-1:1)
    call get_neighbor_cells(flag, ngbr)
    num_neighbor_cells = count(ngbr.eq.1)
  end function num_neighbor_cells

#else

  pure subroutine get_neighbor_cells (flag, ngbr)
    integer, intent(in) :: flag
    integer, intent(inout) :: ngbr(-1:1,-1:1,-1:1)
    integer :: i, j, k
    do k = -1,1
       do j = -1,1
          do i = -1,1
             if (btest(flag,pos_ngbr(i,j,k))) then
                ngbr(i,j,k) = 1
             else
                ngbr(i,j,k) = 0
             end if
          end do
       end do
    end do
  end subroutine get_neighbor_cells

  elemental integer function num_neighbor_cells (flag)
    integer, intent(in) :: flag
    integer ngbr(-1:1,-1:1,-1:1)
    call get_neighbor_cells(flag, ngbr)
    num_neighbor_cells = count(ngbr.eq.1)
  end function num_neighbor_cells

#endif

end module amrex_ebcellflag_module
