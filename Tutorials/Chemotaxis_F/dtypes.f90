module dtypes_module

  use iso_c_binding
  use multifab_module

  implicit none

  integer, parameter :: dim = 2

  integer, parameter :: iu = 1  ! index of population density u
  integer, parameter :: iv = 2  ! index of chemical signal v

  ! chemotaxis context
  type, bind(c) :: cht_ctx_t

     integer(c_int) :: &
          nc = 2, &             ! number of components
          ng = 2, &             ! number of ghost cells
          n_cell = 128          ! number of grid cells

     real(c_double) :: dt, dx, invdx

     real(c_double) :: prob_lo(dim), prob_hi(dim)

     real(c_double) :: &
          diff = 0.1d0, &       ! diffusion constant
          chi  = 5.0d0, &       ! sensitivity constant
          alpha = 0.5d0, &      ! receptor binding constant
          gamma = 10.0d0, &     ! volume filling constant
          phi = 1.0d0           ! saturating chemical production constant
     
  end type cht_ctx_t

  ! options
  type :: cht_opts_t

     character(len=8) :: method = "sdc" ! integration method: "rk" or "sdc"

     integer :: &
          nsteps = 1000, &      ! number of steps
          plot_int = 100        ! plot interval

  end type cht_opts_t

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_namelists(unitno, ctx, opts)
    integer,          intent(in)    :: unitno
    type(cht_ctx_t),  intent(inout) :: ctx
    type(cht_opts_t), intent(inout) :: opts

    character(len=8) :: method
    integer :: n_cell, nsteps, plot_int
    real(8) :: dt

    namelist /options/ method, nsteps, plot_int
    namelist /parameters/ n_cell, dt

    ! read options
    method = opts%method
    nsteps = opts%nsteps
    plot_int = opts%plot_int
    read(unit=unitno, nml=options)
    opts%method = method
    opts%nsteps = nsteps
    opts%plot_int = plot_int

    ! read parameters
    n_cell = ctx%n_cell
    dt = ctx%dt
    read(unit=unitno, nml=parameters)
    ctx%dt = dt
    ctx%n_cell = n_cell

  end subroutine read_namelists

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module dtypes_module
