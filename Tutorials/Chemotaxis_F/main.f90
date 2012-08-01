program main
  
  use boxlib
  use bl_io_module
  use multifab_module
  use chemotaxis_module
  use sdcquad_module


  implicit none

  integer :: un, narg
  character(len=128) :: inputs_file_name

  type(cht_ctx_t)  :: ctx
  type(cht_opts_t) :: opts
  type(sdcquad)    :: sdc

  real(dp_t) :: start_time, run_time, max_run_time

  call boxlib_initialize()
  start_time = parallel_wtime()

  ! read input files
  narg = command_argument_count()
  if ( narg >= 1 ) then
     call get_command_argument(1, value=inputs_file_name)
     un   = unit_new()
     open(unit=un, file=inputs_file_name, status='old', action='read')
     call read_namelists(un, ctx, opts)
     if (opts%method == "sdc") then
        call build(sdc, un)
        call mk_imex_smats(sdc)
     end if
     close(unit=un)
  end if

  ! run chemotaxis solver
  call chemotaxis(ctx, opts, sdc)

  ! deallocate temporary boxarrays and communication mappings, display memory stats
  call layout_flush_copyassoc_cache()

  if ( parallel_ioprocessor() ) then
     print*, " "
     print*, "MEMORY STATS AT END OF PROGRAM"
     print*, " "
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

  ! collect run_time from each processor and store the maximum
  run_time = parallel_wtime() - start_time

  call parallel_reduce(max_run_time, run_time, MPI_MAX, &
                       proc = parallel_ioprocessornode())

  if ( parallel_ioprocessor() ) then
     print*, " "
     print*, "Run time (s) =", max_run_time
  end if

  call boxlib_finalize()

end program main
