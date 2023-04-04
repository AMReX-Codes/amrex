if ( NOT (AMReX_SPACEDIM EQUAL 3) )
   return ()
endif ()

set(_sources     main.cpp)
set(_input_files inputs  particle_file.init  fixed_grids.init)

setup_test(_sources _input_files)

unset(_sources)
unset(_input_files)
