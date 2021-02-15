if (NOT CMAKE_Fortran_COMPILER_LOADED OR NOT (AMReX_SPACEDIM EQUAL 3))
   return()
endif ()

set(_sources MyKernel_F.F90 MyKernel_F.H main.cpp )
set(_inpute_files)

setup_tutorial(_sources _input_files HAS_FORTRAN_MODULE)

unset(_sources)
unset(_input_files)
