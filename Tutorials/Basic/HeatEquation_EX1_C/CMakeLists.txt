if (DIM EQUAL 1)
   return()
endif ()

get_filename_component ( DIRNAME ${CMAKE_CURRENT_LIST_DIR} NAME )

set ( TESTNAME ${DIRNAME} )
set ( EXENAME  ${DIRNAME}.exe )
set ( SRC_DIR  ${CMAKE_CURRENT_LIST_DIR}/Source )
set ( EXE_DIR  ${CMAKE_CURRENT_LIST_DIR}/Exec )

#
# Create target (executable) for this tutorial
# 
add_executable ( ${EXENAME} EXCLUDE_FROM_ALL "" )

target_sources ( ${EXENAME}
   PRIVATE
   ${SRC_DIR}/advance_${DIM}d.f90
   ${SRC_DIR}/advance.cpp
   ${SRC_DIR}/init_phi_${DIM}d.f90
   ${SRC_DIR}/main.cpp
   ${SRC_DIR}/myfunc_F.H
   ${SRC_DIR}/myfunc.H )

set_target_properties ( ${EXENAME} PROPERTIES      
   INCLUDE_DIRECTORIES
   "${SRC_DIR};${CMAKE_CURRENT_BINARY_DIR}/mod_files"
   Fortran_MODULE_DIRECTORY
   ${CMAKE_CURRENT_BINARY_DIR}/mod_files
   RUNTIME_OUTPUT_DIRECTORY
   ${CMAKE_CURRENT_BINARY_DIR} )

target_link_libraries ( ${EXENAME} amrex )


#
# Copy input files to corresponding build dir
# First copy operation is run when builddir is setup
set (INPUT_FILES ${EXE_DIR}/inputs_${DIM}d )

file ( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

add_custom_target ( run_${EXENAME}
   COMMAND
   cp ${INPUT_FILES} ${CMAKE_CURRENT_BINARY_DIR}
   COMMAND
   ${EXENAME} inputs_${DIM}d  nsteps = 10
   DEPENDS
   ${INPUT_FILES}
   WORKING_DIRECTORY
   ${CMAKE_CURRENT_BINARY_DIR}
   COMMENT
   "Running ${EXENAME}"
   VERBATIM
   )


#
# Add to the "tutorial" target
# 
#add_tutorial (${EXENAME})  
add_tutorial (run_${EXENAME})




