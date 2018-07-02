cmake_minimum_required (VERSION 3.5)


project (AMReX)

enable_language (C)
enable_language (CXX)
enable_language (Fortran)


#
# Load required modules 
#
set ( AMREX_CMAKE_MODULES_PATH "${CMAKE_CURRENT_LIST_DIR}/Tools/CMake" CACHE INTERNAL "" )
set ( CMAKE_MODULE_PATH ${AMREX_CMAKE_MODULES_PATH} )

include ( AMReX_Utils )

include ( AMReX_Options )

include ( AMReX_Machines )

#
# Set CMAKE_<LANG>_FLAGS_<CONFIG> if not already defined
# 
set_default_config_flags ()

#
# Provide a default install directory
#
if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
   set ( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/installdir"
      CACHE PATH "AMReX installation directory" FORCE)
endif ()

# 
# Set variable for AMReX versioning
#
find_package (Git QUIET)

set ( TMP "" )

if ( EXISTS ${CMAKE_SOURCE_DIR}/.git AND ${GIT_FOUND} )
   execute_process ( COMMAND git describe --abbrev=12 --dirty --always --tags
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      OUTPUT_VARIABLE TMP )
   string ( STRIP ${TMP} TMP )
endif ()

set ( AMREX_GIT_VERSION "${TMP}" CACHE INTERNAL "" )
unset (TMP)

# 
# Source files for all binaries and libraries found under src
# 
add_subdirectory (Src)

#
# Tutorials
#
add_subdirectory (Tutorials)
