#[=======================================================================[
AMReXBuildInfo
----------------

Provides the function "generate_buildinfo(_target _git_dir)".

generate_buildinfo(_target _git_dir) adds AMReXBuildInfo.H and
AMReXBuildInfo.cpp to the source list of ${_target}.

AMReXBuildInfo.H is a header file located in /path/to/amrex/Tools/C_Scripts
and contains the signatures of functions designed to retrieve build info
like compiler used, source directory paths, compilation flags and so on.

AMReXBuildInfo.cpp is a source file containing the implementations of
the prototypes defined in AMReXBuildInfo.H.
AMReXBuildInfo.cpp is created from AMReXBuildInfo.cpp.in by
generate_buildinfo(_target _git_dir) via the CMake function
configure_file().
AMReXBuildInfo.cpp.in is located in /path/to/amrex/Tools/CMake.


This module

* sets internal variables AMREX_C_SCRIPTS_DIR, AMREX_BUILDINFO_IFILE,
  and AMREX_TOP_DIR

* defines cache variable AMREX_BUILD_DATETIME

* provides the function "generate_buildinfo(_target _git_dir)"

#]=======================================================================]

#
# This include brings in function "get_target_properties_flattened()"
#
include(AMReXTargetHelpers)

#
# Set paths
#
string(REPLACE "/Tools/CMake" "" AMREX_TOP_DIR ${CMAKE_CURRENT_LIST_DIR})
set( AMREX_TOP_DIR ${AMREX_TOP_DIR} CACHE INTERNAL "Top level AMReX directory")

set( AMREX_BUILDINFO_IFILE ${CMAKE_CURRENT_LIST_DIR}/AMReX_buildInfo.cpp.in
   CACHE INTERNAL "Full path and name of AMReX_buildInfo.cpp.in")

set( AMREX_C_SCRIPTS_DIR "${AMREX_TOP_DIR}/Tools/C_scripts"
   CACHE INTERNAL "Path to AMReX' C_scripts dir")

set(AMREX_BUILD_DATETIME "" CACHE STRING
   "User defined build date and time. Set ONLY for reproducibly built binary distributions")

#
# FUNCTION: generate_buildinfo(_target _git_dir)
#
# Adds  AMReXBuildInfo.H and AMReXBuildInfo.cpp to the source
# list of ${_target}.
#
# If ${_git_dir} is a path to a valid git directory, pass
# this info to the build info script for GIT version retrieval.
#
function (generate_buildinfo _target _git_dir)

   #
   # check if _target is a valid target
   #
   if (NOT TARGET ${_target} )
      message(FATAL_ERROR "Target ${_target} does not exists")
   endif ()

   #
   # Set variables to be used in generated source file
   #

   # Date and time -- C++ preprocessor macros
   # If AMREX_BUILD_DATETIME is not set, use macros __DATE__ and __TIME__
   if (AMREX_BUILD_DATETIME)
      set(BUILD_DATE "\"${AMREX_BUILD_DATETIME}\"")
   else()
      set(BUILD_DATE "__DATE__ \" \"  __TIME__" )
   endif ()

   # Target build directory
   get_target_property(BUILD_DIR ${_target} BINARY_DIR)

   # Build machine
   cmake_host_system_information( RESULT BUILD_MACHINE
      QUERY OS_NAME HOSTNAME OS_RELEASE OS_VERSION OS_PLATFORM )
   string(REPLACE ";" " " BUILD_MACHINE "${BUILD_MACHINE}")


   #
   # Build flags
   #
   get_target_properties_flattened(${_target} _includes _defines _flags _link_line)
   get_property(_lang GLOBAL PROPERTY ENABLED_LANGUAGES)
   set(_prop _includes _defines _flags _link_line )

   # Loop over all combinations of language and property and extract
   # what we need
   foreach( _l IN LISTS _lang )

      foreach( _p IN LISTS _prop )

         # _${_l}${_p} is a variable named as _lang_property,
         set(_${_l}${_p} "${${_p}}")
         eval_genex( _${_l}${_p} ${_l} ${CMAKE_${_l}_COMPILER_ID}
           COMP_VERSION ${CMAKE_${_l}_COMPILER_VERSION}
           CONFIG       ${CMAKE_BUILD_TYPE}
           INTERFACE    BUILD)

         if (_${_l}${_p})
            list(REMOVE_DUPLICATES _${_l}${_p})

            if ("${_p}" STREQUAL "_defines")
               string(REPLACE ";" " -D" _${_l}${_p} "-D${_${_l}${_p}}")
            elseif ("${_p}" STREQUAL "_includes")
               string(REPLACE ";" " -I" _${_l}${_p} "-I${_${_l}${_p}}")
            else()
               string(REPLACE ";" " " _${_l}${_p} "${_${_l}${_p}}")
            endif ()
         endif ()

         # Manually escape quotes
         string(REPLACE "\"" "\\\""  _${_l}${_p} "${_${_l}${_p}}")
      endforeach()

      # The configure file needs ${_l}_FLAGS: we set these variables here
      set(${_l}_FLAGS "${CMAKE_${_l}_FLAGS} ${_${_l}_flags}")
      if (CMAKE_BUILD_TYPE)
         string(TOUPPER ${CMAKE_BUILD_TYPE} _ubuild_type)
         set(${_l}_FLAGS "${CMAKE_${_l}_FLAGS_${_ubuild_type}} ${${_l}_FLAGS}")
      endif ()
      set(${_l}_FLAGS "${_${_l}_defines} ${_${_l}_includes} ${${_l}_FLAGS}")
      string(STRIP  "${${_l}_FLAGS}" ${_l}_FLAGS)

   endforeach ()

   # Link flags (we leave it empty and just provide the full link
   # line later
   set(LINK_FLAGS "")

   # Provide full link line instead of just libraries
   # Use c++ link line
   set(LIBRARIES "${_CXX_link_line}")

   # Number of modules -- set to zero for now
   set(NUM_MODULES "int num_modules = 0;")

   #
   # Git hashes for both app code and AMReX if available
   #
   find_package(Git)

   # App code hash
   if (Git_FOUND AND (EXISTS ${_git_dir}/.git))
      # Check the hash of the app code first
      execute_process(
         COMMAND git describe --abbrev=12 --dirty --always --tags
         WORKING_DIRECTORY ${_git_dir}
         OUTPUT_VARIABLE _hash )
      string(STRIP "${_hash}" _hash)
      set(GIT_DECLS "static const char HASH1[] = \"${_hash}\";\n")
      set(GIT_CASE  "case 1: return HASH1;\n")
   endif ()

   # AMReX app code
   if (AMReX_GIT_VERSION)  # If using AMReX as a library
      set(GIT_DECLS "${GIT_DECLS}  static const char HASH2[] = ${AMReX_GIT_VERSION};")
      set(GIT_CASE  "${GIT_CASE}    case 2: return HASH2;")
   elseif (AMREX_GIT_VERSION)  # If using AMReX via add_subdirectory()
      set(GIT_DECLS "${GIT_DECLS}  static const char HASH2[] = \"${AMREX_GIT_VERSION}\";")
      set(GIT_CASE  "${GIT_CASE}    case 2: return HASH2;")
   elseif ( Git_FOUND AND (EXISTS ${AMREX_TOP_DIR}/.git) )   # Backup case
      execute_process(
         COMMAND git describe --abbrev=12 --dirty --always --tags
         WORKING_DIRECTORY ${AMREX_TOP_DIR}
         OUTPUT_VARIABLE _hash )
      string(STRIP "${_hash}" _hash)
      set(GIT_DECLS "${GIT_DECLS}  static const char HASH2[] = \"${_hash}\";")
      set(GIT_CASE  "${GIT_CASE}    case 2: return HASH2;")
   endif ()

   # Other variables to be left unfilled for now
   set(BUILDGIT_DECLS "static const char HASH[] = \"\";")
   set(BUILDGIT_NAME  "static const char NAME[] = \"\";")

   # Generate AMReX_buildInfo.cpp
   configure_file( ${AMREX_BUILDINFO_IFILE}
      ${CMAKE_CURRENT_BINARY_DIR}/AMReX_buildInfo.cpp @ONLY)

   target_sources( ${_target}
      PRIVATE
      ${CMAKE_CURRENT_BINARY_DIR}/AMReX_buildInfo.cpp
      )

   target_sources( ${_target}
      PRIVATE
      ${AMREX_C_SCRIPTS_DIR}/AMReX_buildInfo.H
      )

   target_include_directories( ${_target}
      PUBLIC
      $<BUILD_INTERFACE:${AMREX_C_SCRIPTS_DIR}>
      )

   # Make AMReX_buildInfo.cpp always out of date before building
   add_custom_command( TARGET ${_target}
      PRE_BUILD
      COMMAND ${CMAKE_COMMAND} -E touch_nocreate AMReX_buildInfo.cpp
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )

endfunction ()
