if (NOT CMAKE_Fortran_COMPILER_LOADED OR NOT AMReX_LINEAR_SOLVERS)
   return()
endif ()

set(_sources main.cpp ff.f90)
set(_input_files inputs)

setup_tutorial(_sources _input_files)

unset(_sources)
unset(_input_files)
