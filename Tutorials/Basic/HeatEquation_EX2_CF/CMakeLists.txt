#
# This tutorial needs Fortran
#
if ( NOT CMAKE_Fortran_COMPILER_LOADED )
   return()
endif ()

# List of source files
set(_sources main.cpp advance.cpp myfunc_F.H myfunc.H init_phi_${DIM}d.f90 advance_${DIM}d.f90)
list(TRANSFORM _sources PREPEND "Source/")

# List of input files
file( GLOB_RECURSE _input_files LIST_DIRECTORIES false Exec/input* )

setup_tutorial(_sources _input_files)

unset( _sources )
unset( _input_files   )
