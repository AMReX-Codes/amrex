cmake_minimum_required (VERSION 3.5)


project (AMReX)

# 
# Define the languages used by the project
# 
enable_language (C)
enable_language (CXX)
enable_language (Fortran)


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
# Source files for all binaries and libraries found under src
# 
add_subdirectory (Src)

#
# Tutorials
#
add_subdirectory (Tutorials)
