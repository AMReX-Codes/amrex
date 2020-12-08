#
#
# FUNCTION: generate_config_header
#
# Generate AMReX_Config.H by using the defines set in
# AMREX_DEFINES property
#
#
function ( generate_config_header )

   # Extract from amrex the content of AMREX_DEFINES property
   # This property accumulates the macros that must be defined
   # in the config file
   get_target_property(_defines_list amrex AMREX_DEFINES)
   eval_genex(_defines_list NONE NONE INTERFACE BUILD)

   # set variables from list of defines
   # configure_file() will use these variables to
   # define or undefine items in AMReX_Config.H
   foreach(_define IN LISTS _defines_list)
      string(FIND "${_define}" "=" idx)
      if (idx EQUAL -1) # define without value
         set(${_define} " ")
      else ()   # define with value
         string(SUBSTRING "${_define}" 0      ${idx}  _def)
         string(SUBSTRING "${_define}" ${idx} -1      _val)
         string(REPLACE "=" " " _val "${_val}")
         string(STRIP "${_val}" _val)
         set(${_def} ${_val})
      endif ()

   endforeach ()

   set(COMP_DECLS)
   if (NOT AMReX_DIFFERENT_COMPILER)
       if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
           set(COMPILER_ID_MACRO __GNUC__)
       elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
           set(COMPILER_ID_MACRO __INTEL_COMPILER)
       elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray" )
           set(COMPILER_ID_MACRO  __CRAYC)
       elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI" )
           set(COMPILER_ID_MACRO  __PGI)
       elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
           set(COMPILER_ID_MACRO  __llvm__)
       elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" )
           set(COMPILER_ID_MACRO  __llvm__)
       elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC" )
           set(COMPILER_ID_MACRO  _MSC_VER)
       else ()
           message(FATAL_ERROR "Compiler '${CMAKE_CXX_COMPILER_ID}' not supported by AMReX developers! "
               "Try to configure with -DAMReX_DIFFERENT_COMPILER=ON")
       endif ()
       set(msg "libamrex was built with ${CMAKE_CXX_COMPILER_ID}. To avoid this error, reconfigure with -DAMReX_DIFFERENT_COMPILER=ON")
       set(COMP_DECLS "\n#ifndef ${COMPILER_ID_MACRO}\nstatic_assert(false,\"${msg}\");\n#endif")
   endif()

   configure_file(${AMREX_CMAKE_MODULES_PATH}/AMReX_Config.H.in
      "${AMReX_BINARY_DIR}/AMReX_Config.H")

   target_sources(amrex PRIVATE ${AMReX_BINARY_DIR}/AMReX_Config.H)
   target_include_directories(amrex
      PUBLIC
      $<BUILD_INTERFACE:${AMReX_BINARY_DIR}>)

   install(FILES "${AMReX_BINARY_DIR}/AMReX_Config.H" DESTINATION include)

endfunction ()



#
#
# FUNCTION: add_amrex_define
#
# Populate custom property AMREX_DEFINES with compile definitions.
#
# Arguments:
#
#    _new_define   = variable containing the definition to add.
#                    The new define can be any string with no "-D" prepended.
#                    If the new define is in the form AMREX_SOMENAME,
#                    this function also adds BL_SOMENAME to the list,
#                    unless NO_LEGACY is specified (see below)
#
#    NO_LEGACY     = if specified, the legacy version of a new_define given in the
#                    form AMREX_SOMENAME will not be added.
#
#    IF <cond-var> = new_define is added only if <cond-var> is true
#
#
function ( add_amrex_define _new_define )

   cmake_parse_arguments( "" "NO_LEGACY" "IF" ""  ${ARGN} )

   # Return if flags does not need to be included
   if (_IF AND NOT ${_IF})
      return()
   endif ()

   set(_new_defines ${_new_define})

   # Add legacy definition if required
   if ( NOT _NO_LEGACY )
      string( FIND ${_new_define} "AMREX_" _out )
      if (${_out} GREATER -1 )
	 string(REPLACE "AMREX_" "BL_" _legacy_define ${_new_define})
         list(APPEND _new_defines ${_legacy_define})
      endif ()
   endif ()

   # Populate AMREX_DEFINES (custom) property for amrex target
   get_target_property(_current_defines amrex AMREX_DEFINES)
   if (_current_defines)
      set_target_properties(amrex PROPERTIES AMREX_DEFINES "${_current_defines};${_new_defines}")
   else ()
      set_target_properties(amrex PROPERTIES AMREX_DEFINES "${_new_defines}")
   endif ()

endfunction ()
