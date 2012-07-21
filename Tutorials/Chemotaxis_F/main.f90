program main
  
  use boxlib
  use multifab_module
  use chemotaxis_module
  use bl_io_module

  implicit none

  integer :: un, farg, narg
  character(len=128) :: inputs_file_name

  type(cht_params_t) :: params
  type(cht_opts_t)   :: options

  real(dp_t) :: start_time, run_time, max_run_time

  call boxlib_initialize()
  start_time = parallel_wtime()

  ! read input files
  narg = command_argument_count()
  farg = 1
  if ( narg >= 1 ) then
     call get_command_argument(farg, value=inputs_file_name)
     farg = farg + 1
     un   = unit_new()
     open(unit=un, file=inputs_file_name, status='old', action='read')
     call cht_read_inputs(un, params, options)
     close(unit=un)
  end if

  ! run chemotaxis solver
  call cht_main(params, options)

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
