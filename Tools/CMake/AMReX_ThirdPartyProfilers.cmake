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
function (set_amrex_profilers)

   # 
   # Check if target "amrex" has been defined before
   # calling this macro
   #
   if (NOT TARGET amrex)
      message (FATAL_ERROR "Target 'amrex' must be defined before calling function 'set_amrex_profilers'" )     
   endif ()

   #
   # Check that needed options have already been defined
   # 
   if ( ( NOT ( DEFINED TP_PROFILE ) ) OR ( NOT (DEFINED SITE) ) )
      message ( AUTHOR_WARNING "Required options are not defined" )
      return ()
   endif ()


   if ( ${TP_PROFILE} MATCHES "CRAYPAT" )

      add_amrex_define( AMREX_CRAYPAT )
      if ( ${SITE} MATCHES "nersc" )
         target_include_directories ( amrex PUBLIC $ENV{CRAYPAT_ROOT}/include )
      endif ()

   elseif ( ${TP_PROFILE} MATCHES "FORGE" )

      add_amrex_define( AMREX_FORGE )
      if ( ${SITE} MATCHES "nersc" )
         target_include_directories ( amrex PUBLIC 
	    $ENV{ALLINEA_TOOLS_DIR}/$ENV{ALLINEA_TOOLS_VERSION}/map/wrapper )
      endif()

   elseif ( ${TP_PROFILE} MATCHES "VTUNE" )

      add_amrex_define( AMREX_VTUNE )
      target_compile_options ( amrex PUBLIC -debug inline-debug-info -parallel-source-info=2 )

      if ( ${SITE} MATCHES "nersc" )
         target_compile_options ( amrex PUBLIC -dynamic )
         target_include_directories ( amrex PUBLIC $ENV{VTUNE_AMPLIFIER_XE_2018_DIR}/include )
         target_link_libraries ( amrex PUBLIC -L$ENV{VTUNE_AMPLIFIER_XE_2018_DIR}/lib64 -littnotify )
      endif ()

   endif ()
   
endfunction ()




