module box_util_module

  use box_module
  use boxarray_module
  use ml_boxarray_module

  implicit none

contains

  subroutine read_a_mglib_grid_ba(ba, bx, str)
    use bl_IO_module
    use bl_error_module
    type(boxarray), intent(out) :: ba
    character(len=*), intent(in) :: str
    type(box), intent(out) :: bx
    integer :: un, n
    type(box), allocatable :: bxs(:)
    integer :: i
 
    un = unit_new()
    open(unit=un, file=str, status = 'old', action = 'read', err = 100)
    go to 101
100 call bl_error("READ_A_MGLIB_GRID_BA: failed to open: ", str)
101 continue
    call box_read(bx, un)
    read(unit=un, fmt=*) n
    allocate(bxs(n))
    do i = 1, n
       call box_read(bxs(i), un)
    end do
    call boxarray_build_v(ba, bxs)
    close(unit=un)

  end subroutine read_a_mglib_grid_ba

  subroutine read_a_mglib_grid(mba, str)
    use bl_IO_module
    type(ml_boxarray), intent(out) :: mba
    character(len=*), intent(in) :: str

    call build(mba, 1)
    call read_a_mglib_grid_ba(mba%bas(1), mba%pd(1), str)
    call ml_boxarray_alloc_rr(mba, get_dim(mba%bas(1)))

  end subroutine read_a_mglib_grid

  subroutine read_a_hgproj_grid(mba, str, max_lev_of_mba)
    use bl_IO_module
    use bl_error_module
    type(ml_boxarray), intent(out) :: mba
    character(len=*), intent(in) :: str
    integer, optional, intent(in) :: max_lev_of_mba
    integer :: un, n, nl
    type(box) :: bx1
    type(box), allocatable:: bxs(:)
    integer :: i, j
    integer :: tml
    tml = Huge(tml); if ( present(max_lev_of_mba) ) tml = max_lev_of_mba
    un = unit_new()
    open(unit=un, file=str, status = 'old', action = 'read', err = 100)
    go to 101
100 call bl_error("READ_A_HGPROJ_GRID: failed to open: ", str)
101 continue
    read(unit=un, fmt=*) nl
    nl = min(tml, nl)
    call build(mba, nl)
    do i = 1, nl
       call box_read(bx1, un)
       if ( i == 1 ) then
          call ml_boxarray_alloc_rr(mba, bx1%dim)
          if (nl > 1) mba%rr = 2
       end if
       mba%pd(i) = bx1
       if (i > 1) then
          do n = 1, mba%dim
             mba%rr(i-1,n) = extent(mba%pd(i),n) / extent(mba%pd(i-1),n)
          end do
       end if
       read(unit=un, fmt=*) n
       allocate(bxs(n))
       do j = 1, n
          call box_read(bxs(j), un)
       end do
       call boxarray_build_v(mba%bas(i), bxs)
       deallocate(bxs)
    end do
    close(unit=un)

  end subroutine read_a_hgproj_grid

  subroutine write_a_hgproj_grid(mba, str)
    use bl_IO_module
    use bl_error_module
    type(ml_boxarray), intent(in) :: mba
    character(len=*), intent(in) :: str
    integer :: un
    integer :: i, j
    un = unit_new()
    open(unit=un, file=str, status = 'unknown', action = 'write', err = 100)
    go to 101
100 call bl_error("WRITE_A_HGPROJ_GRID: failed to open: ", str)
101 continue

    write(unit=un,fmt=*) mba%nlevel
    do i = 1, mba%nlevel
       call box_print(mba%pd(i), legacy=.true., advance="NO", unit=un)
       write(unit=un,fmt=*) nboxes(mba, i)
       do j = 1, nboxes(mba, i)
          call box_print(get_box(mba, i, j), legacy=.true., unit=un)        
       enddo
    enddo
    close(unit=un)

  end subroutine write_a_hgproj_grid


  !! Reads a multilevel boxarray in the following format
  !! Record 1: dimension
  !! Record 2: Number of Levels
  !! Record 3: If NL > 1 ; Refinement ratio array(NL-1) values
  !! Record 4: Number of boxes, nb_1, in the coarsest level
  !! Record 5: Box 1 in 1st level lo(1:dm) hi(1:dm)
  !!           repeat nb_1 times.
  !! Goto   4: until number of levels is complete
  !! I/O errors will make this routine fail, as well as allocation
  !! errors, no checks are done.
  subroutine ml_boxarray_read_boxes(mba, str)
    use bl_IO_module
    use bl_error_module
    type(ml_boxarray), intent(out) :: mba
    character(len=*), intent(in) :: str
    integer :: nl, nr, dm, un
    integer, allocatable :: rr(:)
    integer :: i, j
    type(box), allocatable :: bxs(:)
    un = unit_new()
    open(unit=un, file=str, status='old', action = 'read', err = 100)
    go to 101
100 call bl_error("READ_A_MGLIB_GRID_BA: failed to open: ", str)
101 continue
    read(unit=un, fmt=*) dm
    read(unit=un, fmt=*) nl
    call build(mba, nl, dm)
    allocate(rr(nl-1))
    if ( nl > 1 ) then
       read(unit=un, fmt=*) rr
       do i = 1, nl - 1
          mba%rr(i,:) = rr(i)
       end do
    end if
    do i = 1, nl
       read(unit=un, fmt=*) nr
       allocate(bxs(nr))
       do j = 1, nr
          bxs(j)%dim = dm
          read(unit=un, fmt=*) bxs(j)%lo(1:dm), bxs(j)%hi(1:dm)
       end do
       call boxarray_build_v(mba%bas(i), bxs)
       deallocate(bxs)
    end do
    mba%pd(1) = boxarray_bbox(mba%bas(1))
    do i = 2, nl
       mba%pd(i) = refine(mba%pd(i-1), mba%rr(i-1,:))
    end do
    close(unit=un)
  end subroutine ml_boxarray_read_boxes

end module box_util_module
