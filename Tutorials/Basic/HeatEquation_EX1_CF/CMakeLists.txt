if ( (DIM EQUAL 1) OR NOT ENABLE_FORTRAN)
   return()
endif ()

set(_sources
   advance_${DIM}d.f90
   advance.cpp
   init_phi_${DIM}d.f90
   main.cpp
   myfunc_F.H
   myfunc.H )

list(TRANSFORM _sources PREPEND "${CMAKE_CURRENT_LIST_DIR}/Source/")

set(_include_dirs ${CMAKE_CURRENT_LIST_DIR}/Source)
set( _input_files ${CMAKE_CURRENT_LIST_DIR}/Exec/inputs_${DIM}d)

setup_tutorial(_sources _include_dirs _input_files HAS_FORTRAN_MODULES)

unset( _sources )
unset( _include_dirs )
unset( _input_files  )



# if ( (NOT ENABLE_CUDA) AND (NOT ENABLE_DPCPP) )

#    add_custom_target( run_${EXENAME}
#       COMMAND
#       cp ${INPUT_FILES} ${CMAKE_CURRENT_BINARY_DIR}
#       COMMAND
#       ${EXENAME} inputs_${DIM}d  nsteps = 10
#       DEPENDS
#       ${INPUT_FILES}
#       WORKING_DIRECTORY
#       ${CMAKE_CURRENT_BINARY_DIR}
#       COMMENT
#       "Running ${EXENAME}"
#       VERBATIM
#       )



#    #
#    # Add to the "tutorial" target
#    #
#    add_tutorial (run_${EXENAME})

# endif ()
