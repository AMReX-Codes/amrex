cmake_minimum_required (VERSION 3.5)


project (AMReX)

# 
# Define the languages used by the project
# 
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
set ( CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/Tools/CMake )

#
# Include module with custom macros and functions
#
include ( AMReX_Utils )

#
# Initialize all the global variables
#
include ( AMReX_CMakeVariables )

#
# Load AMReX options 
#
include ( AMReX_Options )

#
# Configure build
# 
include ( AMReX_Config )

# 
# Source files for all binaries and libraries found under src
# 
add_subdirectory (Src)

#
# Tutorials
#
add_subdirectory (Tutorials)
