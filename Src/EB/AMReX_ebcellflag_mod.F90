
module amrex_ebcellflag_module
  use amrex_fort_module, only : amrex_real
  implicit none
  private
  public :: is_regular_cell, is_single_valued_cell, is_multi_valued_cell, &
       is_covered_cell, get_neighbor_cells, num_neighbor_cells, &
       set_regular_cell, amrex_ebcellflag_count, get_neighbor_cells_int_single
  
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

#endif

  elemental subroutine set_regular_cell (flag)
    integer, intent(inout) :: flag
    call mvbits(regular, 0, w_type, flag, 0)
  end subroutine set_regular_cell

  subroutine amrex_ebcellflag_count(lo,hi,flag,flo,fhi,nregular,nsingle,nmulti,ncovered) &
       bind(c,name='amrex_ebcellflag_count')
    integer, dimension(3), intent(in) :: lo(3), hi(3), flo(3), fhi(3)
    integer, intent(in) :: flag(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    integer, intent(out) :: nregular,nsingle,nmulti,ncovered
    integer :: i,j,k
    nregular = 0
    nsingle = 0
    nmulti = 0
    ncovered = 0
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (is_regular_cell(flag(i,j,k))) then
                nregular = nregular+1
             else if (is_single_valued_cell(flag(i,j,k))) then
                nsingle = nsingle+1
             else if (is_multi_valued_cell(flag(i,j,k))) then
                nmulti = nmulti+1
             else
                ncovered = ncovered+1
             end if
          end do
       end do
    end do
  end subroutine amrex_ebcellflag_count

end module amrex_ebcellflag_module
