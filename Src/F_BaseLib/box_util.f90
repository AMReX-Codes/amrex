module box_util_module

  use box_module
  use boxarray_module
  use mboxarray_module

  implicit none

contains

  !! Places a box into the world.
  !! The center of the is uniformly distributed in the world.
  !! The size of the box is uniformly distrubed between MN and MX
  !! The resulting box is intersected with the world
  function box_random_box(world, mn, mx, mt) result(r)
    use mt19937_module
    type(mt19937), intent(inout), optional :: mt
    type(box) :: r
    type(box), intent(in) :: world
    integer, intent(in) :: mn, mx
    real(kind=dp_t) :: aa(world%dim,2)
    integer ::  spot(world%dim)
    integer :: hwide(world%dim)

    if ( present(mt) ) then
       call mt_random_number(mt, aa)
    else
       call mt_random_number(aa)
    end if
    hwide = (mn + aa(:,1)*(mx-mn))/2
    spot = lwb(world) + aa(:,2)*(upb(world)-lwb(world))
    call build(r, spot-hwide, spot+hwide)
    r = intersection(r, world)

  end function box_random_box

  !! Makes an array of random boxes
  subroutine make_random_boxes_bv(bxs, world, mn, mx, mt)
    use mt19937_module
    type(box), dimension(:), intent(out) :: bxs
    type(mt19937), intent(inout), optional :: mt
    integer, intent(in) :: mn, mx
    type(box), intent(in) :: world
    integer i

    do i = 1, size(bxs)
       bxs(i) = box_random_box(world, mn, mx, mt)
    end do

  end subroutine make_random_boxes_bv

  subroutine read_a_mglib_grid_ba(ba, bx, str)
    use bl_IO_module
    type(boxarray), intent(out) :: ba
    character(len=*), intent(in) :: str
    type(box), intent(out) :: bx
    integer :: un, n
    type(box), allocatable :: bxs(:)
    integer :: i
 
    un = unit_new()
    open(unit=un, file=str, status = 'old', action = 'read')
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
    type(mboxarray), intent(out) :: mba
    character(len=*), intent(in) :: str

    call build(mba, 1)
    call read_a_mglib_grid_ba(mba%bas(1), mba%pd(1), str)
    call mboxarray_alloc_rr(mba, mba%bas(1)%dim)

  end subroutine read_a_mglib_grid

  subroutine read_a_hgproj_grid(mba, str)
    use bl_IO_module
    type(mboxarray), intent(out) :: mba
    character(len=*), intent(in) :: str
    integer :: un, n, nl
    type(box) :: bx1
    type(box), allocatable:: bxs(:)
    integer :: i, j
 
    un = unit_new()
    open(unit=un, file=str, status = 'old', action = 'read')
    read(unit=un, fmt=*) nl
    call build(mba, nl)
    do i = 1, nl
       call box_read(bx1, un)
       if ( i == 1 ) then
          call mboxarray_alloc_rr(mba, bx1%dim)
          mba%rr = 2
       end if
       mba%pd(i) = bx1
       if (i > 1) mba%rr(i-1,:) = box_extent_d(mba%pd(i),1) / box_extent_d(mba%pd(i-1),1)
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
  subroutine mboxarray_read_boxes(mba, str)
    use bl_IO_module
    type(mboxarray), intent(out) :: mba
    character(len=*), intent(in) :: str
    integer :: nl, nr, dm, un
    integer, allocatable :: rr(:)
    integer :: i, j
    type(box), allocatable :: bxs(:)
    un = unit_new()
    open(unit=un, file=str, status='old', action = 'read')
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
  end subroutine mboxarray_read_boxes

  subroutine build_random_boxarray(ba, pd, nb, mn, mx, bf)
    use mt19937_module
    implicit none
    type(boxarray), intent(out) :: ba
    integer, intent(in) :: nb, mn, mx, bf
    type(box), intent(in) :: pd
    type(box) :: tpd
    type(box) :: bxs(nb)

    integer :: tmn, tmx

    tpd = coarsen(pd, bf)
    tmn = max(mn/bf, 1)
    tmx = max(mx/bf, 1)

    call make_random_boxes_bv(bxs, tpd, tmn, tmx)
    call boxarray_build_v(ba, bxs)
    call boxarray_refine(ba, bf)

    call boxarray_to_domain(ba)

  end subroutine build_random_boxarray

end module box_util_module
