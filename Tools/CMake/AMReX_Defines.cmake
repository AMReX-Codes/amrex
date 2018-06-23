#
#
#  This file provides:
#
#    AMREX_DEFINES:          list of the cpp flags to be used on ALL sources
#    AMREX_Fortran_DEFINES:  list of Fortran-specific cpp flags
#    add_define():           function to add definitions to AMREX_DEFINES 
#
#  Once this file is included, AMREX_DEFINES and AMREX_Fortran_DEFINES will be
#  populated with the preprocessor directives needed for a succesfull build.
#  Further CPP flags can be appended manually or via add_define() after this
#  file has been included.
# 
#  As per xSDK requirements, if a user set the env variable CPPFLAGS,
#  CPPFLAGS should overwrite AMReX_DEFINES. Since this is not possible without
#  breaking the build (most of the defines here are necessary for AMReX to compile),
#  for the time being we will not follow this requirement.
#  
# 

#
# Function to accumulate preprocessor directives (to be removed)
#
function ( add_define new_define )

   cmake_parse_arguments ( DEFINE "" "IF" ""  ${ARGN} )

   set ( condition  1 )

   if (DEFINE_IF)
      set ( condition ${${DEFINE_IF}} )
   endif()

   # Return if flags does not need to be included
   if (NOT condition)
      return ()
   endif ()
   
   set ( definition -D${new_define} )
   
   list ( FIND AMREX_DEFINES ${definition} out )
   
   if ( ${out} EQUAL -1 )
      list ( APPEND AMREX_DEFINES ${definition} )
      #  Add legacy definition
      string (FIND ${definition} "AMREX_" out )
      if (${out} GREATER -1 )
	 string (REPLACE "AMREX_" "BL_" legacy_definition ${definition})
	 list ( APPEND AMREX_DEFINES ${legacy_definition} )
      endif ()
      set ( AMREX_DEFINES ${AMREX_DEFINES} PARENT_SCOPE )
   endif ()
   
endfunction ()


function ( set_amrex_defines )

   # 
   # Check if target "amrex" has been defined before
   # calling this macro
   #
   if (NOT TARGET amrex)
      message (FATAL_ERROR "Target 'amrex' must be defined before calling function 'set_amrex_defines'" )
   endif ()

   #
   # Set compile definition property
   #
   
   # Git version
   add_amrex_define ( "AMREX_GIT_VERSION=\"${AMREX_GIT_VERSION}\"" NO_LEGACY )

   # XSDK mode 
   add_amrex_define ( AMREX_XSDK IF USE_XSDK_DEFAULTS )

   # Debug flag
   add_amrex_define ( AMREX_DEBUG NO_LEGACY IF DEBUG )
   if ( (NOT DEBUG) AND (NOT USE_XSDK_DEFAULTS) )
      add_amrex_define (NDEBUG)
   endif ()

   # Base profiling options
   add_amrex_define ( AMREX_PROFILING IF ENABLE_BASE_PROFILE )

   add_amrex_define ( AMREX_TRACE_PROFILING IF ENABLE_TRACE_PROFILE )

   add_amrex_define ( AMREX_COMM_PROFILING  IF ENABLE_COMM_PROFILE )

   # Tiny profiler
   add_amrex_define ( AMREX_TINY_PROFILING IF ENABLE_TINY_PROFILE )

   # Mem profiler 
   add_amrex_define ( AMREX_MEM_PROFILING IF ENABLE_MEM_PROFILE )

   # Profparser 
   add_amrex_define ( AMREX_USE_PROFPARSER IF ENABLE_PROFPARSER )

   # Backtrace
   if (ENABLE_BACKTRACE)
      add_amrex_define ( AMREX_BACKTRACING )
      add_amrex_define ( AMREX_TESTING )
   endif ()

   # Third party profiling
   if (${TP_PROFILE} MATCHES "CRAYPAT")
      add_amrex_define ( AMREX_CRAYPAT )
   elseif (${TP_PROFILE} MATCHES "FORGE")
      add_amrex_define ( AMREX_FORGE )
   elseif (${TP_PROFILE} MATCHES "VTUNE")
      add_amrex_define ( AMREX_VTUNE )
   endif ()

   # MPI
   add_amrex_define ( AMREX_USE_MPI IF ENABLE_MPI )

   # OpenMP
   add_amrex_define ( AMREX_USE_OMP IF ENABLE_OMP )

   # Precision
   if (NOT ENABLE_DP)
      add_amrex_define (AMREX_USE_FLOAT)
   endif ()

   # Fortran/C mangling scheme
   add_amrex_define ( AMREX_FORT_USE_${FORTLINK} )

   # Dimensionality
   add_amrex_define ( AMREX_SPACEDIM=${DIM} )

   # System
   add_amrex_define ( AMREX_${CMAKE_SYSTEM_NAME} )

   # Particles
   if (ENABLE_PARTICLES)
      add_amrex_define ( AMREX_PARTICLES )

      if (NOT ENABLE_DP_PARTICLES)
	 add_amrex_define ( AMREX_SINGLE_PRECISION_PARTICLES )
      endif ()
   endif ()

   #  Assertions
   add_amrex_define ( AMREX_USE_ASSERTION IF ENABLE_ASSERTIONS )

   add_amrex_define ( AMREX_USE_EB IF ENABLE_EB )

   add_amrex_define ( AMREX_USE_F_INTERFACES IF ENABLE_FORTRAN_INTERFACES )
   
   add_amrex_define ( AMREX_NO_STRICT_PREFIX )

   #
   # Fortran-specific defines: BL_LANG_FORT e AMREX_LANG_FORT
   #
   target_compile_definitions ( amrex PUBLIC
      $<$<COMPILE_LANGUAGE:Fortran>:BL_LANG_FORT;AMREX_LANG_FORT> )
   
endfunction ()


function ( add_amrex_define new_define )

   # 
   # Check if target "amrex" has been defined before
   # calling this macro
   #
   if (NOT TARGET amrex)
      message (FATAL_ERROR "Target 'amrex' must be defined before calling function 'add_amrex_define'" )
   endif ()
   
   cmake_parse_arguments ( DEFINE "NO_LEGACY" "IF" ""  ${ARGN} )

   set ( condition  1 )

   if (DEFINE_IF)
      set ( condition ${${DEFINE_IF}} )
   endif()

   # Return if flags does not need to be included
   if (NOT condition)
      return ()
   endif ()

   target_compile_definitions ( amrex PUBLIC ${new_define} )

   if ( NOT DEFINE_NO_LEGACY )
      # Add legacy definition
      string ( FIND ${new_define} "AMREX_" out )

      if (${out} GREATER -1 )
	 string (REPLACE "AMREX_" "BL_" legacy_define ${new_define})
	 target_compile_definitions ( amrex PUBLIC ${legacy_define} )
      endif ()
   endif () 
   
endfunction ()



