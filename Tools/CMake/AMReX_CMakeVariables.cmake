#
# 
#  This file contains all the global variables necessary
#  to the build process. Some of them are only defined
#  here and will be assigned a value someplace else.
#  Some other variables, for example  AMREX_DEFAULT_INSTALL_DIR,
#  are initialized with a specific value which is meant to be used
#  unchanged.
# 
#  

#
# Check if project() has been called; abort if not
# 
if ( NOT PROJECT_NAME )
   message ( FATAL_ERROR "AMReX_CMakeVariables.cmake cannot be included\
before calling project()" )
endif ()

#
# Check if this file has been loaded already
#
if ( DEFINED __AMREX_CMAKEVARIABLES__ )
   return ()
endif ()

set (  __AMREX_CMAKEVARIABLES__ "" )

# Set paths for build system
set ( CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod_files )

# The type of build ( will need to be uppercase )
set ( AMREX_BUILD_TYPE )

# Provide a default install directory
set ( AMREX_DEFAULT_INSTALL_DIR "${PROJECT_SOURCE_DIR}/installdir" )

if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
   set ( CMAKE_INSTALL_PREFIX "${AMREX_DEFAULT_INSTALL_DIR}" CACHE PATH
      "AMReX installation directory" FORCE)
endif ()


# 
# Set variable for AMReX versioning
#
find_package (Git QUIET)
   
# Check whether .git is present and git installed and, if so,
# retrieve infos
set ( output "" )
if ( EXISTS ${CMAKE_SOURCE_DIR}/.git AND ${GIT_FOUND} )
   execute_process ( COMMAND git describe --abbrev=12 --dirty --always --tags
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      OUTPUT_VARIABLE output )
   string ( STRIP ${output} output )
endif ()
set ( AMREX_GIT_VERSION ${output} )

# 
# Set directory paths
# 
set (AMREX_SOURCE_DIR         ${CMAKE_SOURCE_DIR}/Src )
set (AMREX_CMAKE_MODULES_DIR  ${CMAKE_CURRENT_LIST_DIR})
set (AMREX_LIBRARY_DIR        ${CMAKE_INSTALL_PREFIX}/lib )
set (AMREX_INCLUDE_DIR        ${CMAKE_INSTALL_PREFIX}/include )
set (AMREX_BINARY_DIR         ${CMAKE_INSTALL_PREFIX}/bin )
set (AMREX_TOOLS_DIR          ${CMAKE_INSTALL_PREFIX}/Tools )
set (AMREX_CMAKE_DIR          ${CMAKE_INSTALL_PREFIX}/cmake )
set (AMREX_LIBRARIES          amrex )

# 
# Config files for export
# 
set ( AMREX_CONFIG_INSTALL_INFILE  ${AMREX_CMAKE_MODULES_DIR}/AMReXConfig.cmake.in)
set ( AMREX_CONFIG_INSTALL_OUTFILE ${PROJECT_BINARY_DIR}/AMReXConfig.cmake)

#
# Load host system info
#
include ( AMReX_Machines )

# 
# Compiler flags 
# 
set ( AMREX_Fortran_FLAGS )
set ( AMREX_C_FLAGS )
set ( AMREX_CXX_FLAGS )

#
# Compile- and link-time variables
#
#    AMREX_EXTRA_<LANG>_FLAGS :
#               list of flags needed by external packages
#    AMREX_EXTRA_<LANG>_INCLUDE_PATH :
#               list of include paths needed by external packages
#    AMREX_EXTRA_<LANG>_LINK_LINE :
#               link line to append at link time to account for external packages
#
set (AMREX_EXTRA_C_INCLUDE_PATH)
set (AMREX_EXTRA_CXX_INCLUDE_PATH)
set (AMREX_EXTRA_Fortran_INCLUDE_PATH)
set (AMREX_EXTRA_C_FLAGS)
set (AMREX_EXTRA_CXX_FLAGS)
set (AMREX_EXTRA_Fortran_FLAGS)
set (AMREX_EXTRA_C_LIBRARIES)
set (AMREX_EXTRA_CXX_LIBRARIES)
set (AMREX_EXTRA_Fortran_LIBRARIES)
set (AMREX_EXTRA_C_LINK_FLAGS)
set (AMREX_EXTRA_CXX_LINK_FLAGS)
set (AMREX_EXTRA_Fortran_LINK_FLAGS)
set (AMREX_EXTRA_C_LINK_LINE)
set (AMREX_EXTRA_CXX_LINK_LINE)
set (AMREX_EXTRA_Fortran_LINK_LINE)


#
# Detect Fortran name mangling scheme for C/Fortran interface and
# set FORTLINK variable accordingly
#
include ( FortranCInterface )
include ( ${FortranCInterface_BINARY_DIR}/Output.cmake )

set ( FORTLINK "" )

if ( FortranCInterface_GLOBAL_SUFFIX STREQUAL "" )

   set (FORTLINK "${FortranCInterface_GLOBAL_CASE}CASE" )
   message (STATUS "Fortran name mangling scheme: ${FORTLINK} (no append underscore)")

elseif ( (FortranCInterface_GLOBAL_SUFFIX STREQUAL "_")  AND
      ( FortranCInterface_GLOBAL_CASE STREQUAL "LOWER" ) )

   set (FORTLINK "UNDERSCORE")
   message (STATUS "Fortran name mangling scheme: ${FORTLINK} (lower case, append underscore)")

else ()
   
   message (AUTHOR_WARNING "Fortran to C mangling not compatible with AMReX code")
   
endif ()
