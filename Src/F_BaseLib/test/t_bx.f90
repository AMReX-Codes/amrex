subroutine t_bx
  use box_module
  implicit none
  type(box) :: bx
  bx = allbox(2)
  print *, bx
  bx = allbox_grc(bx, grow = 2)
  print *, bx
contains
  function allbox_grc(bx, grow, refine, coarsen) result(r)
    type(box) :: r
    type(box), intent(in) :: bx
    integer, intent(in), optional :: grow, refine, coarsen
    integer :: dm
    r = bx
    dm = r%dim
    if ( present(grow) ) then
    else if ( present(refine) ) then
    else if ( present(coarsen) ) then
    end if
  end function allbox_grc
end subroutine t_bx

subroutine t_ba_self_intersection
  use ml_boxarray_module
  use box_util_module
  use bl_prof_module
  implicit none
  type(boxarray) :: ba
  integer :: dm
  integer :: ng
  integer :: pd_xyz(MAX_SPACEDIM)
  integer :: ba_maxsize
  logical pmask(MAX_SPACEDIM)
  type(ml_boxarray) :: mba
  character(len=64) :: test_set
  integer :: i, f, n
  type(box) :: bx
  type(bl_prof_timer), save :: bpt, bpt_r, bpt_s
  
  call build(bpt, "t_ba_self_intersection")

  ng = 1
  test_set = "grids.5034"

  call build(bpt_r, "ba_read")
  call read_a_mglib_grid(mba, test_set)
  call destroy(bpt_r)

  ba = mba%bas(1)

  call build(bpt_s, "ba_s")
  do i = 1, nboxes(ba)
     bx = grow(get_box(ba,i),1)
     call self_intersection(bx, ba)
  end do
  call destroy(bpt_s)

  call destroy(mba)
  call destroy(bpt)
contains

  subroutine self_intersection(bx, ba)
    type(box), intent(in) :: bx
    type(boxarray), intent(in) :: ba
    integer :: i
    type(bl_prof_timer), save :: bpt_i
    type(box) :: bx1
    call build(bpt_i, "ba_i")
    do i = 1, nboxes(ba)
       bx1 = intersection(bx, get_box(ba,i))
    end do
    call destroy(bpt_i)

  end subroutine self_intersection

end subroutine t_ba_self_intersection

