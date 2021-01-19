if ( (AMReX_SPACEDIM EQUAL 1) OR (NOT AMReX_FORTRAN_INTERFACES) )
    return()
endif ()

set(_sources main.F90 mytest.F90 init_prob.F90)

file( GLOB_RECURSE _input_files LIST_DIRECTORIES false ${CMAKE_CURRENT_LIST_DIR}/input*)

setup_tutorial(_sources _input_files HAS_FORTRAN_MODULES)

unset(_sources)
unset(_input_files)
