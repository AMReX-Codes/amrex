
subroutine fmain () bind(c)

  use baselib_module

  implicit none

  integer :: n_cell, max_grid_size, nsteps, plot_int
  type(ParmParse) :: pp
  type(Box) :: domain
  type(BoxArray) :: bs
  type(Geometry) :: geom

  ! ParmParse is way of reading inputs from the inputs file
  ! "get" means it must be set in the inputs file, whereas
  ! "query" means it may not may not be in the inputs file
  call parmparse_build(pp)

  call pp%get("n_cell", n_cell);  ! # of cells in each dimension
  call pp%get("nsteps", nsteps) ! # of steps

  max_grid_size = 32   ! default max grid size
  call pp%query("max_grid_size", max_grid_size);

  plot_int = -1 ! default to no plotfiles
  call pp%query("plot_int", plot_int);
  
  ! Define a single box covering the domain
  domain = Box((/0,0,0/), (/n_cell-1, n_cell-1, n_cell-1/))

  ! Initialize the boxarray "bs" from the single box "bx"
  call boxarray_build(bs, domain)

  ! Break up boxarray "bs" into chunks no larger than "max_grid_size" along a direction
  call bs%maxSize(max_grid_size)

  ! This defines a Geometry object.
  call geometry_build(geom, domain)


end subroutine fmain
