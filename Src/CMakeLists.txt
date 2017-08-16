#
# Add all the needed subdirs to the include path
#
set ( EXCLUDE_PATHS )  # Path to exclude from include dirs

if (NOT ENABLE_FBASELIB)
   list (APPEND EXCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/F_BaseLib)
   # The includes belove are for Fortran-based MG
   list (APPEND EXCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/LinearSolvers/C_to_F_MG)
   list (APPEND EXCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/LinearSolvers/F_MG)
endif ()

if (NOT ENABLE_LINEAR_SOLVERS)
   list (APPEND EXCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/LinearSolvers)
endif ()

if (NOT ENABLE_FORTRAN_INTERFACES) 
   list (APPEND EXCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/F_Interfaces)
endif ()

if (NOT ENABLE_PARTICLES)
   list (APPEND EXCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/Particle)
endif ()

list (APPEND EXCLUDE_PATHS ${CMAKE_CURRENT_LIST_DIR}/EBAMRElliptic) # not yet included in build 


find_include_paths (INCLUDE_PATHS EXCLUDE ${EXCLUDE_PATHS} ) 


include_directories ( ${CMAKE_Fortran_MODULE_DIRECTORY} ${INCLUDE_PATHS}
   ${AMREX_EXTRA_Fortran_INCLUDE_PATH} ${AMREX_EXTRA_CXX_INCLUDE_PATH} )

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
include (AmrCore/CMakeLists.txt)
include (Amr/CMakeLists.txt)
include (Boundary/CMakeLists.txt)
include (Extern/CMakeLists.txt)

#
# Optional components
#
if (ENABLE_FBASELIB)
   include (F_BaseLib/CMakeLists.txt)
endif ()


if (ENABLE_LINEAR_SOLVERS)
   include (LinearSolvers/CMakeLists.txt)
endif ()

if (ENABLE_FORTRAN_INTERFACES) 
   include (F_Interfaces/CMakeLists.txt)
endif ()

if ( ENABLE_PARTICLES )
   include (Particle/CMakeLists.txt)
endif ()


#
# Define install  
#
set_target_properties ( amrex PROPERTIES PUBLIC_HEADER "${AMREX_ALLHEADERS}" )

install ( TARGETS amrex EXPORT amrex ARCHIVE DESTINATION ${AMREX_LIBRARY_DIR}
   PUBLIC_HEADER DESTINATION ${AMREX_INCLUDE_DIR} )

install ( DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${AMREX_INCLUDE_DIR})


#
# This header in a weird path has to be copied to install includes
# 
install ( FILES ${PROJECT_SOURCE_DIR}/Tools/C_scripts/AMReX_buildInfo.H
   DESTINATION ${AMREX_INCLUDE_DIR} )

#
# Write and install configure file 
# 
configure_file ( ${AMREX_CONFIG_INSTALL_INFILE} ${AMREX_CONFIG_INSTALL_OUTFILE} @ONLY)
install ( FILES ${AMREX_CONFIG_INSTALL_OUTFILE} DESTINATION ${AMREX_CMAKE_DIR})

#
# Install Tools directory
#
install (DIRECTORY ${PROJECT_SOURCE_DIR}/Tools
   DESTINATION ${CMAKE_INSTALL_PREFIX}
   USE_SOURCE_PERMISSIONS )

