cmake_minimum_required (VERSION 3.3)

project (AMReX)

# Define the languages used by the project
enable_language (C)
enable_language (CXX)
enable_language (Fortran)

#
# Require C++11 standard
#
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON) 
set (CMAKE_CXX_EXTENSIONS OFF)

#
# Always verbose output during make
# 
set ( CMAKE_VERBOSE_MAKEFILE on )


#
# Set the search path for cmake modules
#
set ( AMREX_TOOLS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Tools )
set ( CMAKE_MODULE_PATH ${AMREX_TOOLS_DIR}/CMake )

#
# Initialize all the global variables
#
include ( AMReX_CMakeVariables )

#
# Include modules with custom macros and functions
#
include ( AMReX_Utils )

########## Load CMake tools directory ##########
# Useful variables pointing to directories in the source tree
set(CCSE_TOOLS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Tools")
set(CCSE_PERL_DIR "${CCSE_TOOLS_DIR}/C_scripts")
set(CCSE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}")

# Set the module search path so find_package and include commands
# can locate files in <root source tree>/tools/cmake
set(CCSE_MODULE_PATH "${CCSE_TOOLS_DIR}/CMake")

# Additional build options
include ( AMReX_Version )
include ( AMReX_Options )
include ( AMReX_Config )

#
# Usage:
#
# Typically cmake is used to support an "out-of-source-tree" build.  To this with default
# settings (see below):
#
#     #> export bdir=`pwd`
#     #> mkdir <new_build_dir>
#     #> cd <new_build_dir>
#     #> cmake -DCMAKE_INSTALL_PREFIX:PATH=<some_install_prefix_for_amrex> $bdir
#     #> make
#     #> make install
#
# The install target will place all C++/Fortran header and mod files in an include folder,
# the CCSE libraries in a lib folder, plus some stuff in perl and cmake folders.
#
# An example cmake input file to build an app that links to the libs built as above:
#
#              ./Tutorials/MultiGrid_C/CMakeLists.txt
#
# To build this tutorial with cmake, follow instructions similar to the above for building
# AMReX:
#
#     #> mkdir <new_build_dir>
#     #> cd <new_build_dir>
#     #> cmake -DAMREX_INSTALL_PREFIX:PATH=<some_install_prefix_for_amrex> $bdir/Tutorials/MultiGrid_C
#     #> make
#     #> mpiexec -np 4 ./mgc_tutorial inputs.3d
#

# Define the project name
# This command will define
#     CCSE_SOURCE_DIR
#     CCSE_BINARY_DIR


# None, Debug, Release, RelWithDebInfo, MinSizeRel
#set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of CMake build")

set(CMAKE_INSTALL_PREFIX "${CURRENT_CMAKE_BINDARY_DIR}" CACHE STRING "Root dir where AMReX files installed")

# Enable testing, ctest needs this
# all add_test commands are ignored unless this is called!
enable_testing()

# A property for accumulating the a global CCSE link line.
set_property(GLOBAL PROPERTY CCSE_LINK_LINE "-L${CMAKE_INSTALL_PREFIX}/lib")

# A property for accumulating CCSE library targets
set_property(GLOBAL PROPERTY CCSE_LIBRARY_TARGETS)
include(InstallManager)


# Source files for all binaries and libraries found under src
add_subdirectory(Src)

# Source file for postprocess exes
add_subdirectory(Tools/Postprocessing/F_Src) 

# Make the exports only _after_ doing the build
create_exports()

install_include_files (${CCSE_PERL_DIR}/AMReX_buildInfo.H)
ADD_INSTALL_CMAKE_FILES(${CCSE_MODULE_PATH}/CCSEOptions.cmake)
ADD_INSTALL_CMAKE_FILES(${CCSE_MODULE_PATH}/FindCCSE.cmake)
ADD_INSTALL_CMAKE_FILES(${CCSE_MODULE_PATH}/CCSEConfigVersion-install.cmake.in)


option(ENABLE_Config_Report "Print out detailed information at the end of a configuration")
set(CCSE_CONFIG_LOG "${CCSE_BINARY_DIR}/ccse-config.log"
    CACHE string "CCSE configuration log file")
# include(CCSEConfigReport)
# if (ENABLE_Config_Report)
#     execute_process(COMMAND "${cat_exec}" "${CCSE_CONFIG_LOG}" OUTPUT_VARIABLE config_output)
#     message(STATUS "contents of: ${CCSE_CONFIG_LOG}")
#     message(STATUS "${config_output}")
# endif()    
