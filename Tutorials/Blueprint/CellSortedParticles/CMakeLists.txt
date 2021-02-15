if (NOT (AMReX_SPACEDIM EQUAL 3) OR NOT CMAKE_Fortran_COMPILER_LOADED OR NOT AMReX_PARTICLES)
   return()
endif ()

set(_sources cell_sorted_3d.F90 cell_sorted_F.H CellSortedPC.cpp CellSortedPC.H main.cpp)
set(_input_files inputs)

setup_tutorial(_sources _input_files)

unset(_sources)
unset(_input_files)
