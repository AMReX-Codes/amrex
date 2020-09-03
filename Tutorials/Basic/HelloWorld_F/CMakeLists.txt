#
# This tutorial needs Fortran
#
if ( NOT CMAKE_Fortran_COMPILER_LOADED )
   return()
endif ()

# List of source files
set(_sources fmain.f90)

# List of include directories
set(_include_dirs)

# List of input files
set(_input_files)

setup_tutorial(_sources _include_dirs _input_files)

unset( _sources )
unset( _include_dirs  )
unset( _input_files   )
