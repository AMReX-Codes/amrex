if ( NOT (3 IN_LIST AMReX_SPACEDIM) OR AMReX_MPI
     OR AMReX_PRECISION STREQUAL "SINGLE"
     OR AMReX_PARTICLES_PRECISION STREQUAL "SINGLE" )
   return ()
endif ()

set(_sources     main.cpp)
set(_input_files inputs  particle_file.init  fixed_grids.init)

# A smaller problem keeps the ctest short.
setup_test(3 _sources _input_files
   CMDLINE_PARAMS_AFTER_INPUTS nx=64 ny=64 nz=64)

unset(_sources)
unset(_input_files)
