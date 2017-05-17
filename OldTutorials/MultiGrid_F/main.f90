program main

    ! These modules are all in BoxLib/Src/F_BaseLib
    use boxlib
    use multifab_module
    use bl_IO_module

    ! These modules are all in this directory
    use constants
    use init_rhs_mod
    use init_coeffs_mod
    use cycle_tower_mod
    use solve_mod

    implicit none

    ! From inputs file
    integer          :: n_cells
    integer          :: max_boxsize
    integer          :: dim
    integer          :: s1
    integer          :: s2
    integer          :: v_bot
    integer          :: fmg_bot
    integer          :: max_iter
    integer          :: cycle_type
    integer          :: interp_type
    integer          :: verbose        ! verbose > 2 will print arrays, so be careful with large data sets
    logical          :: memory_verbose ! this controls the printing of the memory stats only and is set below
    integer          :: rhs_type
    integer          :: coeffs_type
    double precision :: eps

    ! Other local variables
    integer                       :: cyc, nlvl
    integer, allocatable          :: lo(:), hi(:)
    double precision              :: dx ! Grid spacing
    double precision, allocatable :: prob_lo(:), prob_hi(:) ! Problem domain
    logical, allocatable          :: is_periodic(:) ! Periodic boundary conditions
    type(cycle_tower)             :: ct ! Holds all information for multigrid solver
    type(box)                     :: bx
    type(boxarray)                :: ba
    type(layout)                  :: la
    type(multifab)                :: rhs, coeffs

    ! Dummy variables used to read inputs file
    integer :: un, farg, narg
    logical :: need_inputs_file, found_inputs_file
    character(len=128) :: fname

    ! Used to keep track of program run time
    double precision :: start_time, run_time, run_time_IOproc

    namelist /probin/ n_cells, max_boxsize, dim, s1, s2, v_bot, fmg_bot, max_iter, cycle_type, eps, &
                      interp_type, rhs_type, coeffs_type, verbose

    ! If running in parallel. this will print out the number of MPI processes
    ! and OpenMP threads
    call boxlib_initialize()

    memory_verbose = .false.

    ! parallel_wtime() returns the number of wallclock-time seconds since
    ! the program began
    start_time = parallel_wtime()

    ! Default values -- will be overwritten by inputs file
    n_cells = 64
    max_boxsize = 32
    dim = 2
    s1 = 2
    s2 = 2
    v_bot = 8
    fmg_bot = 40
    max_iter = 50
    cycle_type = 1
    eps = 1E-10
    interp_type = 1
    rhs_type = 1
    coeffs_type = 1
    verbose = 1

    ! Read inputs file and overwrite any default values
    narg = command_argument_count()
    need_inputs_file = .true.
    farg = 1
    if ( need_inputs_file .AND. narg >= 1 ) then
        call get_command_argument(farg, value = fname)
        inquire(file = fname, exist = found_inputs_file )
        if (found_inputs_file) then
            farg = farg + 1
            un = unit_new()
            open(unit=un, file = fname, status = 'old', action = 'read')
            read(unit=un, nml = probin)
            close(unit=un)
            need_inputs_file = .false.
        end if
    end if

    ! Easy way to modify one input at a time
    do while (farg <= narg)
        call get_command_argument(farg, value=fname)
        select case(fname)
            case('--n_cells')
                farg = farg + 1
                call get_command_argument(farg, value=fname)
                read(fname,*) n_cells
            case('--max_boxsize')
                farg = farg + 1
                call get_command_argument(farg, value=fname)
                read(fname,*) max_boxsize
            case('--dim')
                farg = farg + 1
                call get_command_argument(farg, value=fname)
                read(fname,*) dim
            case('--s1')
                farg = farg + 1
                call get_command_argument(farg, value=fname)
                read(fname,*) s1
            case('--s2')
                farg = farg + 1
                call get_command_argument(farg, value=fname)
                read(fname,*) s2
            case('--v_bot')
                farg = farg + 1
                call get_command_argument(farg, value=fname)
                read(fname,*) v_bot
            case('--fmg_bot')
                farg = farg + 1
                call get_command_argument(farg, value=fname)
                read(fname,*) fmg_bot
            case('--max_iter')
                farg = farg + 1
                call get_command_argument(farg, value=fname)
                read(fname,*) max_iter
            case('--cycle_type')
                farg = farg + 1
                call get_command_argument(farg, value=fname)
                read(fname,*) cycle_type
            case('--eps')
                farg = farg + 1
                call get_command_argument(farg, value=fname)
                read(fname,*) eps
            case('--interp_type')
                farg = farg + 1
                call get_command_argument(farg, value=fname)
                read(fname,*) interp_type
            case('--verbose')
                farg = farg + 1
                call get_command_argument(farg, value=fname)
                read(fname,*) verbose
            case('--rhs_type')
                farg = farg + 1
                call get_command_argument(farg, value=fname)
                read(fname,*) rhs_type
            case('--coeffs_type')
                farg = farg + 1
                call get_command_argument(farg, value=fname)
                read(fname,*) coeffs_type
            case default
                print *, 'ERROR - Unknown input', fname
                stop
        end select
        farg = farg + 1
    end do

    ! Now that we have dim, we can allocate these
    allocate(lo(dim))
    allocate(hi(dim))
    allocate(is_periodic(dim))
    allocate(prob_lo(dim))
    allocate(prob_hi(dim))

    ! Physical problem is from (0,0) to (1,1), periodic on all sides
    prob_lo(:) = 0
    prob_hi(:) = 1.d0
    is_periodic(:) = .true.
    dx = (prob_hi(1)-prob_lo(1))/n_cells

    ! Create a box from (0,0) to (n_cells-1, n_cells-1)
    lo(:) = 0
    hi(:) = n_cells-1
    bx = make_box(lo,hi)

    ! Initialize the boxarray to be a single box
    call boxarray_build_bx(ba,bx)

    ! Overwrite the boxarray to respect max_boxsize
    call boxarray_maxsize(ba,max_boxsize)

    ! Build layout.  The third argument, bx, is the problem domain.
    call layout_build_ba(la,ba,bx,pmask=is_periodic)
    call destroy(ba)

    ! Build multifab for RHS and coefficients using same layout.
    call multifab_build(rhs,la,nc=1,ng=0)
    call multifab_build(coeffs,la,nc=1,ng=1)

    ! Fill them in
    call init_rhs(rhs,rhs_type,prob_lo,dx)
    call init_coeffs(coeffs,coeffs_type)

    ! Number of multigrid levels (each box coarsens down to 2x2)
    nlvl = nint(log(dble(max_boxsize))/log(2.d0))

    ! Build tower object
    call build_cycle_tower(ct,rhs,coeffs,la,dim,n_cells,dx,s1,s2,v_bot,fmg_bot,nlvl,max_iter,&
                           eps,cycle_type,interp_type,rhs_type,coeffs_type,verbose)

    if (parallel_IOProcessor()) then
        print *, 'MULTIGRID SOLVE'
    end if

    ! Solve using multigrid
    call solve(ct)

    ! Make sure to destroy everything or else you could leak memory
    call destroy_cycle_tower(ct)
    call destroy(la)
    deallocate(lo)
    deallocate(hi)
    deallocate(is_periodic)
    deallocate(prob_lo)
    deallocate(prob_hi)

    ! Deallocate temporary boxarrays and communication mappings
    call layout_flush_copyassoc_cache()

    ! Check for memory that should have been deallocated
    if (memory_verbose) then
        if (parallel_IOProcessor()) then
            print *, 'MEMORY STATS AT THE END OF PROGRAM'
            print *, ''
        end if
        call print(multifab_mem_stats(),    "    multifab")
        call print(fab_mem_stats(),         "         fab")
        call print(boxarray_mem_stats(),    "    boxarray")
        call print(layout_mem_stats(),      "      layout")
        call print(boxassoc_mem_stats(),    "    boxassoc")
        call print(fgassoc_mem_stats(),     "     fgassoc")
        call print(syncassoc_mem_stats(),   "   syncassoc")
        call print(copyassoc_mem_stats(),   "   copyassoc")
        call print(fluxassoc_mem_stats(),   "   fluxassoc")
    end if

    ! parallel_wtime() returns the number of wallclock-time seconds since the program began
    run_time = parallel_wtime() - start_time

    ! collect run_time from each processor and store the maximum
    call parallel_reduce(run_time_IOproc, run_time, MPI_MAX, proc=parallel_IOProcessorNode())

    if (parallel_IOprocessor()) then
        print *, 'Run time (sec) = ', run_time_IOproc
    end if

    call boxlib_finalize()

end program main


