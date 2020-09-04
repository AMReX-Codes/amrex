if ( NOT (DIM EQUAL 3) )
   return ()
endif ()

set(_sources cell_sorted_3d.F90 cell_sorted_F.H CellSortedPC.cpp CellSortedPC.H main.cpp )
set(_input_files inputs)

setup_tutorial(_sources _input_files HAS_FORTRAN_MODULES)

unset(_sources)
unset(_input_files)
