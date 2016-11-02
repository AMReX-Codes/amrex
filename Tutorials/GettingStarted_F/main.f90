program main

  use boxlib
  use multifab_module
  use work_on_data_module

  implicit none

  ! box, boxarray, layout, and multifab
  type(box) :: bx
  type(boxarray) :: ba
  type(layout) :: la
  type(multifab) :: data

  ! dimensionality of problem
  integer :: dm

  ! stuff we will allocate later once we have dm
  integer, allocatable :: lo(:), hi(:)
  logical, allocatable :: is_periodic(:)

  ! if running in parallel, this will print out the number of MPI 
  ! processes and OpenMP threads
  call boxlib_initialize()

  ! here we hard-code this problem to two dimensions
  ! for real applications, we would read this in from an inputs file
  dm = 2

  allocate(lo(dm),hi(dm))
  allocate(is_periodic(dm))

  lo(1:dm) = 0
  hi(1:dm) = 15
  is_periodic(1:dm) = .true.

  bx = make_box(lo,hi)
  call boxarray_build_bx(ba,bx) ! ba has one 16^2 box
  call boxarray_maxsize(ba,8)   ! ba now has four 8^2 boxes

  ! build a boxarray
  call layout_build_ba(la,ba,bx,is_periodic)

  ! build a multifab with 2 components and 6 ghost cells
  call multifab_build(data,la,2,6)

  ! look in work_on_data.f90
  call work_on_data(data)

  ! fill periodic and interior ghost cells
  call multifab_fill_boundary(data)

  ! free up memory to prevent leaks
  call multifab_destroy(data)

  call boxlib_finalize()

end program main
