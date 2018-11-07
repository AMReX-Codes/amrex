cmake_minimum_required (VERSION 3.5)

project (AMReX)

enable_language (C)
enable_language (CXX)
enable_language (Fortran)

#
# Load required modules 
#
set( AMREX_CMAKE_MODULES_PATH "${CMAKE_CURRENT_LIST_DIR}/Tools/CMake" CACHE INTERNAL "" )
set( CMAKE_MODULE_PATH ${AMREX_CMAKE_MODULES_PATH} )

#
# Provide a default install directory
#
if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
   set ( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/installdir"
      CACHE PATH "AMReX installation directory" FORCE)
endif ()

#
# Include options, utilities and other stuff we need
# 
include( AMReX_Utils )
include( AMReX_Options )
include( AMReX_Machines )

#
# Set CMAKE_<LANG>_FLAGS_<CONFIG> if not already defined
# 
set_default_config_flags ()

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
else ()
   # Grep first line from file CHANGES if AMREX_GIT_VERSION is not
   # provided
   file(STRINGS ${PROJECT_SOURCE_DIR}/CHANGES ALL_VERSIONS REGEX "#")
   list(GET ALL_VERSIONS 0 TMP)
   string(REPLACE "#" "" TMP "${TMP}")
   string (STRIP "${TMP}" TMP )
   set (TMP "${TMP}.0")
endif ()

set ( AMREX_GIT_VERSION "${TMP}" CACHE INTERNAL "" )
unset (TMP)

#
# Install Blitz and Algoim at configuration time
# if needed.
# This is only required if ENABLE_3D_NODAL_MLMG is ON
# and external libraries paths are not given 
#
if ( ENABLE_3D_NODAL_MLMG AND
      NOT ( ALGOIM_INSTALL_DIR AND BLITZ_INSTALL_DIR ) )
   include(AMReX_InstallExternalLibs) 
endif ()


# 
# Source files for all binaries and libraries found under src
# 
add_subdirectory (Src)

#
# Tutorials and "test_install" target)
#
add_subdirectory (Tutorials)



