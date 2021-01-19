if ( NOT (AMReX_SPACEDIM EQUAL 3) OR NOT CMAKE_Fortran_COMPILER_LOADED OR AMReX_CUDA OR AMReX_HIP)
   return ()
endif ()

set(_sources cell_sorted_3d.F90 cell_sorted_F.H CellSortedPC.cpp CellSortedPC.H main.cpp )
set(_input_files inputs)

setup_tutorial(_sources _input_files HAS_FORTRAN_MODULES)

unset(_sources)
unset(_input_files)
