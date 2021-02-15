if ( NOT AMReX_FORTRAN_INTERFACES )
   return()
endif ()

# List of source files
set(_sources fmain.f90)

# List of input files
set(_input_files)

setup_tutorial(_sources _input_files)

unset( _sources )
unset( _input_files   )
