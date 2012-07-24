module dtypes_module

  use multifab_module

  integer, parameter :: iu = 1  ! index of population density u
  integer, parameter :: iv = 2  ! index of chemical signal v

  ! chemotaxis context
  type :: cht_ctx_t

     integer :: dim = 2         ! number of dimensions
     integer :: nc = 2          ! number of components
     integer :: ng = 2          ! number of ghost cells
     integer :: n_cell = 32     ! number of grid cells

     type(layout) :: la

     double precision :: dx
     
  end type cht_ctx_t

  ! options
  type :: cht_opts_t

     character(len=8) :: method = "sdc" ! integration method: "rk" or "sdc"

  end type cht_opts_t

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cht_read_inputs(unitno, ctx, opts)
    integer,          intent(in)    :: unitno
    type(cht_ctx_t),  intent(inout) :: ctx
    type(cht_opts_t), intent(inout) :: opts

    character(len=8) :: method
    integer :: n_cell

    namelist /options/ method
    namelist /parameters/ n_cell

    ! read options
    method = opts%method
    read(unit=unitno, nml=options)
    opts%method = method

    ! read parameters
    n_cell = ctx%n_cell
    read(unit=unitno, nml=parameters)
    ctx%n_cell = n_cell

  end subroutine cht_read_inputs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module dtypes_module
