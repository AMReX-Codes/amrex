subroutine t_boxassoc_1
  use layout_module

  implicit none
  type(boxassoc) :: bxasc
  type(boxarray) :: ba
  integer :: dm
  integer :: ng
  integer :: pd_xyz(MAX_SPACEDIM)
  integer :: ba_maxsize
  type(layout) :: la
  logical :: nodal(MAX_SPACEDIM)
  type(box) :: bxs(4)
  logical pmask(MAX_SPACEDIM)
  type(box) :: pd
  integer :: i

  dm = 2
  ng = 1
  pd_xyz     = 32
  ba_maxsize = 16

  nodal = .true.
  nodal = .false.
  pmask = .true.

  if ( .false. ) then
     bxs(1) = make_box((/0,0/),(/3,3/))
     bxs(2) = make_box((/4,4/),(/7,7/))
     bxs(3) = make_box((/0,8/),(/7,16/))
     call build(ba, bxs(1:2))
  else
     call build(ba, make_box((/(0,i=1,dm)/), pd_xyz(1:dm)-1))
     call boxarray_maxsize(ba, ba_maxsize)
  end if

  pd = bbox(ba)
  call build(la, ba, pd = pd, pmask = pmask(1:dm))

  call boxassoc_build(bxasc, la%lap, ng, nodal = nodal(1:dm), cross = .false.)

  call boxassoc_destroy(bxasc)
  call destroy(la)
  call boxarray_destroy(ba)

end subroutine t_boxassoc_1

subroutine t_boxassoc
  use layout_module
  use ml_boxarray_module
  use box_util_module
  implicit none
  type(boxassoc) :: bxasc
  type(boxarray) :: ba
  integer :: dm
  integer :: ng
  type(layout) :: la
  logical pmask(MAX_SPACEDIM)
  type(box) :: pd
  type(ml_boxarray) :: mba
  character(len=64) :: test_set
  logical :: nodal(MAX_SPACEDIM)

  ng = 1
  test_set = "grids.5034"

  nodal = .false.

  call read_a_mglib_grid(mba, test_set)

  ba = mba%bas(1)
  pd = mba%pd(1)

  dm = get_dim(ba)

  call build(la, ba, pd = pd, pmask = pmask(1:dm))

  call boxassoc_build(bxasc, la%lap, ng, nodal = nodal, cross = .false.)

  !! call boxassoc_print(bxasc)

  call boxassoc_destroy(bxasc)

  call destroy(la)
  call destroy(mba)

end subroutine t_boxassoc

