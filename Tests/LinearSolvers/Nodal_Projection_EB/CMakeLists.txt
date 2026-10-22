if ( (NOT AMReX_EB) OR (NOT AMReX_LINEAR_SOLVERS_INCFLO) OR NOT (3 IN_LIST AMReX_SPACEDIM))
   return()
endif ()

set(_sources main.cpp)
set(_input_files inputs_3d)

# A smaller problem keeps the ctest short.
setup_test(3 _sources _input_files
   CMDLINE_PARAMS_AFTER_INPUTS n_cell=32 max_grid_size=16)

unset(_sources)
unset(_input_files)
