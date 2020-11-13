if (NOT AMReX_FORTRAN_INTERFACES OR NOT AMReX_MPI)
   return()
endif ()

set(_sources main.F90)
set(_input_files)

setup_tutorial(_sources _input_files)

unset(_sources)
unset(_input_files)
