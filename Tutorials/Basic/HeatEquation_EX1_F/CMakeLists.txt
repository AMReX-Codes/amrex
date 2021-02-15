#
# This tutorial needs Fortran interfaces
# Disable for cuda
#
if ( NOT AMReX_FORTRAN_INTERFACES OR AMReX_CUDA )
   return()
endif ()

# List of source files
set(_sources fmain.F90 init_phi.f90 advance.f90)

# List of input files
file( GLOB_RECURSE _input_files LIST_DIRECTORIES false input* )

setup_tutorial(_sources _input_files HAS_FORTRAN_MODULES)

unset( _sources )
unset( _input_files   )
