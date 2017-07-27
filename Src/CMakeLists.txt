#
# Add all the subdir to the include path
# For now exclude the following directory because of the conflicting name between
# header files
#
find_include_paths (INCLUDE_PATHS EXCLUDE ${CMAKE_CURRENT_LIST_DIR}/EBAMRElliptic)

include_directories ( ${CMAKE_Fortran_MODULE_DIRECTORY} ${INCLUDE_PATHS}
   ${AMREX_EXTRA_Fortran_INCLUDE_PATH} ${AMREX_EXTRA_CXX_INCLUDE_PATH} )

install (DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)

#
# Define a macro to accumulate all headers
#
set (AMREX_ALLHEADERS "" )
macro (accumulate_headers headers)
   list (APPEND AMREX_ALLHEADERS ${${headers}} )
endmacro ()

#
# Fist, define the library we want to add
#
add_library ( amrex "")


#
# Now, one by one, let's add all the source files in the subdirectories
#
include (Base/CMakeLists.txt)
include (F_BaseLib/CMakeLists.txt)
include (AmrCore/CMakeLists.txt)
include (Amr/CMakeLists.txt)
include (LinearSolvers/CMakeLists.txt)
include (Boundary/CMakeLists.txt)
include (F_Interfaces/CMakeLists.txt) 
include (Extern/CMakeLists.txt)

if ( ENABLE_PARTICLES )
   include (Particle/CMakeLists.txt)
endif ()

#
# Define the install dir
#
set_target_properties ( amrex PROPERTIES PUBLIC_HEADER "${AMREX_ALLHEADERS}" )
install ( TARGETS amrex
   ARCHIVE DESTINATION ${AMREX_LIBRARY_DIR}
   PUBLIC_HEADER DESTINATION ${AMREX_INCLUDE_DIR}
   )

# Move Tools dir to install location
install(CODE "MESSAGE(\"-- Copying   : ${AMREX_TOOLS_DIR} --> ${CCSE_INSTALL_DIR}\")" )
install(CODE "FILE( COPY ${AMREX_TOOLS_DIR} DESTINATION ${CCSE_INSTALL_DIR} )")
# SET_TARGET_PROPERTIES(cboxlib PROPERTIES LINKER_LANGUAGE C)
# SET_TARGET_PROPERTIES(cfboxlib PROPERTIES LINKER_LANGUAGE C)
# SET_TARGET_PROPERTIES(fboxlib PROPERTIES LINKER_LANGUAGE C)

# add_subdirectory(${CAMRDATADIR})


