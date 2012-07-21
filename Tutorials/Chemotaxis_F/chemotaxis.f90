module chemotaxis_module

  use boxlib
  use multifab_module
  use write_plotfile_module

  implicit none

  ! parameters
  type :: cht_params_t

     integer :: dim = 2         ! number of dimensions
     integer :: nc = 2          ! number of components
     integer :: ng = 2          ! number of ghost cells
     integer :: n_cell = 32     ! number of grid cells
     
  end type cht_params_t

  ! options
  type :: cht_opts_t

     character(len=8) :: method = "sdc" ! integration method: "rk" or "sdc"

  end type cht_opts_t

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cht_read_inputs(unitno, params, opts)
    integer,            intent(in)    :: unitno
    type(cht_params_t), intent(inout) :: params
    type(cht_opts_t),   intent(inout) :: opts

    character(len=8) :: method
    integer :: n_cell

    namelist /options/ method
    namelist /parameters/ n_cell

    ! read options
    method = opts%method
    read(unit=unitno, nml=options)
    opts%method = method

    ! read parameters
    n_cell = params%n_cell
    read(unit=unitno, nml=parameters)
    params%n_cell = n_cell

  end subroutine cht_read_inputs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cht_main(params, opts)
    type(cht_params_t), intent(inout) :: params
    type(cht_opts_t),   intent(inout) :: opts

    integer,    allocatable :: lo(:), hi(:)
    real(dp_t), allocatable :: prob_lo(:), prob_hi(:)
    logical,    allocatable :: is_periodic(:)
    
    real(dp_t)     :: dx, dt, time
  
    type(box)      :: bx
    type(boxarray) :: ba
    type(layout)   :: la
    type(multifab) :: U

    ! allocate
    allocate(lo(params%dim),hi(params%dim))
    allocate(prob_lo(params%dim),prob_hi(params%dim))
    allocate(is_periodic(params%dim))

    ! physical problem is a box on (-1,-1) to (1,1), periodic on all sides
    prob_lo     = -1.d0
    prob_hi     =  1.d0
    is_periodic = .true.

    ! create a box from (0,0) to (n_cell-1,n_cell-1)
    lo = 0
    hi = params%n_cell-1
    bx = make_box(lo,hi)

    dx = (prob_hi(1)-prob_lo(1)) / params%n_cell

    ! build layout
    call build(ba,bx)
    ! call boxarray_maxsize(ba,max_grid_size)
    call build(la,ba,bx,pmask=is_periodic)

    ! build multifab
    call build(U,la,params%nc,params%ng)
    
    ! set initial condition and write plot 0
    call cht_initial(U)
    call write_plotfile(la,U,0,dx,0.0d0,prob_lo,prob_hi)

    ! run
    select case(opts%method)
    case("rk")
       print*, "RK"
    case("sdc")
       print*, "SDC"
    end select

    ! destroy/deallocate
    call destroy(ba)
    call destroy(la)
    call destroy(U)

    deallocate(lo,hi,prob_lo,prob_hi,is_periodic)

  end subroutine cht_main

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cht_initial(U)
    type(multifab), intent(inout) :: U

    
  end subroutine cht_initial

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module chemotaxis_module
