#
# This file provides the following variables
#
# FUNCTION: set_amrex_profilers
#
# Setup the target "amrex" to use a third party profiler
# Before using this function, target "amrex" must have been constructed.
# This function returns right away if global variables TP_PROFILE or SITE
# have not been defined before the call
#
# Author: Michele Rosso
# Date  : June 26, 2018
#
#
function (set_amrex_profilers AMREX_TARGET)

   #
   # Check if target "amrex" has been defined before
   # calling this macro
   #
   if (NOT TARGET ${AMREX_TARGET})
      message (FATAL_ERROR "Target '${AMREX_TARGET}' must be defined before calling function 'set_amrex_profilers'" )
   endif ()

   cmake_host_system_information( RESULT _machine QUERY HOSTNAME )

   if ( TP_PROFILE MATCHES "CRAYPAT" )

      add_amrex_define( AMREX_CRAYPAT )
      if ( _machine MATCHES "cori" )
         target_include_directories(${AMREX_TARGET} PUBLIC $ENV{CRAYPAT_ROOT}/include )
      endif ()

   elseif ( TP_PROFILE MATCHES "FORGE" )

      add_amrex_define( AMREX_FORGE )
      if ( _machine MATCHES "cori" )
         target_include_directories(${AMREX_TARGET} PUBLIC
	    $ENV{ALLINEA_TOOLS_DIR}/$ENV{ALLINEA_TOOLS_VERSION}/map/wrapper )
      endif()

   elseif ( TP_PROFILE MATCHES "VTUNE" )

      add_amrex_define( AMREX_VTUNE )
      target_compile_options(${AMREX_TARGET} PUBLIC -debug inline-debug-info -parallel-source-info=2 )

      if ( _machine MATCHES "cori" )
         target_compile_options(${AMREX_TARGET} PUBLIC -dynamic )
         target_include_directories(${AMREX_TARGET} PUBLIC $ENV{VTUNE_AMPLIFIER_XE_2018_DIR}/include )
         target_link_libraries(${AMREX_TARGET} PUBLIC -L$ENV{VTUNE_AMPLIFIER_XE_2018_DIR}/lib64 -littnotify )
      endif ()

   endif ()

endfunction ()
