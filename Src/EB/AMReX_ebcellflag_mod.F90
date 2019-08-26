
module amrex_ebcellflag_module
  use amrex_fort_module, only : amrex_real
  implicit none
  private
  public :: is_regular_cell, is_single_valued_cell, is_multi_valued_cell, &
       is_covered_cell, get_cell_type, get_neighbor_cells, num_neighbor_cells, pos_ngbr, &
       set_regular_cell, set_covered_cell, set_single_valued_cell, set_neighbor, clear_neighbor, &
       clear_allneighbors, is_neighbor, get_neighbor_cells_int_single, &
       regular, single_valued, multi_valued, covered
  
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

  interface get_neighbor_cells
     module procedure get_neighbor_cells_int
     module procedure get_neighbor_cells_real
  end interface get_neighbor_cells

contains

  elemental function get_neighbor_cells_int_single (flag,i,j,k)
    integer, intent(in) :: flag
    integer, intent(in) :: i, j, k
    integer :: get_neighbor_cells_int_single

    if (btest(flag,pos_ngbr(i,j,k))) then
      get_neighbor_cells_int_single = 1
    else
      get_neighbor_cells_int_single = 0
    end if
  end function get_neighbor_cells_int_single
  
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

  elemental function get_cell_type (flag) result(r)
    integer, intent(in) :: flag
    integer :: r
    r = ibits(flag,0,w_type)
  end function get_cell_type

  elemental function clear_allneighbors (flag) result(r)
    integer, intent(in) :: flag
    integer :: r
    r = flag
    call mvbits(0, 0, 27, r, pos_ngbr(-1,-1,-1))
    r = ibset(r, pos_ngbr(0,0,0))
  end function clear_allneighbors

#if (AMREX_SPACEDIM == 2)

  pure subroutine get_neighbor_cells_int (flag, ngbr)
    integer, intent(in) :: flag
    integer, intent(out) :: ngbr(-1:1,-1:1)
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
  end subroutine get_neighbor_cells_int

  pure subroutine get_neighbor_cells_real (flag, ngbr)
    integer, intent(in) :: flag
    real(amrex_real), intent(out) :: ngbr(-1:1,-1:1)
    integer :: i, j
    do j = -1,1
       do i = -1,1
          if (btest(flag,pos_ngbr(i,j,0))) then
             ngbr(i,j) = 1._amrex_real
          else
             ngbr(i,j) = 0._amrex_real
          end if
       end do
    end do
  end subroutine get_neighbor_cells_real

  elemental integer function num_neighbor_cells (flag)
    integer, intent(in) :: flag
    integer ngbr(-1:1,-1:1)
    call get_neighbor_cells_int(flag, ngbr)
    num_neighbor_cells = count(ngbr.eq.1)
  end function num_neighbor_cells

  elemental function set_neighbor (flag,i,j) result(r)
    integer, intent(in) :: flag,i,j
    integer :: r
    r = ibset(flag, pos_ngbr(i,j,0))
  end function set_neighbor

  elemental function clear_neighbor (flag,i,j) result(r)
    integer, intent(in) :: flag,i,j
    integer :: r
    r = ibclr(flag, pos_ngbr(i,j,0))
  end function clear_neighbor

  elemental function is_neighbor (flag,i,j) result(r)
    integer, intent(in) :: flag,i,j
    logical :: r
    r = btest(flag,pos_ngbr(i,j,0))
  end function is_neighbor

#else

  pure subroutine get_neighbor_cells_int (flag, ngbr)
    integer, intent(in) :: flag
    integer, intent(out) :: ngbr(-1:1,-1:1,-1:1)
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
  end subroutine get_neighbor_cells_int

  pure subroutine get_neighbor_cells_real (flag, ngbr)
    integer, intent(in) :: flag
    real(amrex_real), intent(out) :: ngbr(-1:1,-1:1,-1:1)
    integer :: i, j, k
    do k = -1,1
       do j = -1,1
          do i = -1,1
             if (btest(flag,pos_ngbr(i,j,k))) then
                ngbr(i,j,k) = 1._amrex_real
             else
                ngbr(i,j,k) = 0._amrex_real
             end if
          end do
       end do
    end do
  end subroutine get_neighbor_cells_real

  elemental integer function num_neighbor_cells (flag)
    integer, intent(in) :: flag
    integer ngbr(-1:1,-1:1,-1:1)
    call get_neighbor_cells_int(flag, ngbr)
    num_neighbor_cells = count(ngbr.eq.1)
  end function num_neighbor_cells

  elemental function set_neighbor (flag,i,j,k) result(r)
    integer, intent(in) :: flag,i,j,k
    integer :: r
    r = ibset(flag, pos_ngbr(i,j,k))
  end function set_neighbor

  elemental function clear_neighbor (flag,i,j,k) result(r)
    integer, intent(in) :: flag,i,j,k
    integer :: r
    r = ibclr(flag, pos_ngbr(i,j,k))
  end function clear_neighbor

  elemental function is_neighbor (flag,i,j,k) result(r)
    integer, intent(in) :: flag,i,j,k
    logical :: r
    r = btest(flag,pos_ngbr(i,j,k))
  end function is_neighbor

#endif

  elemental subroutine set_regular_cell (flag)
    integer, intent(inout) :: flag
    call mvbits(regular, 0, w_type, flag, 0)
  end subroutine set_regular_cell

  elemental subroutine set_covered_cell (flag)
    integer, intent(inout) :: flag
    call mvbits(covered, 0, w_type, flag, 0)
    call mvbits(0, 0, w_numvofs, flag, pos_numvofs)
  end subroutine set_covered_cell

  elemental subroutine set_single_valued_cell (flag)
    integer, intent(inout) :: flag
    call mvbits(single_valued, 0, w_type, flag, 0)
  end subroutine set_single_valued_cell

end module amrex_ebcellflag_module
