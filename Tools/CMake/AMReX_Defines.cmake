# 
# FUNCTION: set_amrex_defines
#
# Populate the "compile definitions" property of target "amrex".
# Target "amrex" must exist before calling this function.
#
# As per xSDK requirements, if a user set the env variable CPPFLAGS,
# CPPFLAGS should overwrite AMReX_DEFINES. Since this is not possible without
# breaking the build (most of the defines here are necessary for AMReX to compile),
# for the time being we will not follow this requirement.
#  
# Author: Michele Rosso
# Date  : June 26, 2018
#
# 
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
   add_amrex_define( "AMREX_GIT_VERSION=\"${AMREX_GIT_VERSION}\"" NO_LEGACY )

   # XSDK mode 
   add_amrex_define( AMREX_XSDK IF USE_XSDK_DEFAULTS )

   # Debug flag
   add_amrex_define( AMREX_DEBUG NO_LEGACY IF DEBUG )

   # Base profiling options
   add_amrex_define( AMREX_PROFILING IF ENABLE_BASE_PROFILE )
   add_amrex_define( AMREX_TRACE_PROFILING IF ENABLE_TRACE_PROFILE )
   add_amrex_define( AMREX_COMM_PROFILING  IF ENABLE_COMM_PROFILE )

   # Tiny profiler
   add_amrex_define( AMREX_TINY_PROFILING IF ENABLE_TINY_PROFILE )

   # Mem profiler 
   add_amrex_define( AMREX_MEM_PROFILING IF ENABLE_MEM_PROFILE )

   # Profparser 
   add_amrex_define( AMREX_USE_PROFPARSER IF ENABLE_PROFPARSER )

   # Backtrace
   if (ENABLE_BACKTRACE)
      add_amrex_define( AMREX_BACKTRACING )
      add_amrex_define( AMREX_TESTING )
   endif ()

   # Third party profiling
   if (${TP_PROFILE} MATCHES "CRAYPAT")
      add_amrex_define( AMREX_CRAYPAT )
   elseif (${TP_PROFILE} MATCHES "FORGE")
      add_amrex_define( AMREX_FORGE )
   elseif (${TP_PROFILE} MATCHES "VTUNE")
      add_amrex_define( AMREX_VTUNE )
   endif ()

   # MPI
   add_amrex_define( AMREX_USE_MPI IF ENABLE_MPI )

   # OpenMP
   add_amrex_define( AMREX_USE_OMP IF ENABLE_OMP )

   # CUDA
   add_amrex_define( AMREX_USE_CUDA    IF ENABLE_CUDA )
   add_amrex_define( AMREX_USE_GPU     IF ENABLE_CUDA )
   
   # Precision
   if (NOT ENABLE_DP)
      add_amrex_define(AMREX_USE_FLOAT)
   endif ()

   # Dimensionality
   add_amrex_define( AMREX_SPACEDIM=${DIM} )

   # System
   add_amrex_define( AMREX_${CMAKE_SYSTEM_NAME} )

   # Particles
   if (ENABLE_PARTICLES)
      add_amrex_define( AMREX_PARTICLES )
      if (NOT ENABLE_DP_PARTICLES)
	 add_amrex_define( AMREX_SINGLE_PRECISION_PARTICLES )
      endif ()
   endif ()

   # External libraries for nodal MLMG
   add_amrex_define( USE_ALGOIM IF ENABLE_3D_NODAL_MLMG)
   
   #  Assertions
   add_amrex_define( AMREX_USE_ASSERTION IF ENABLE_ASSERTIONS )
   add_amrex_define( AMREX_USE_EB IF ENABLE_EB )
   add_amrex_define( AMREX_USE_F_INTERFACES IF ENABLE_FORTRAN_INTERFACES ) 
   add_amrex_define( AMREX_NO_STRICT_PREFIX )

   #
   # Fortran-specific defines: BL_LANG_FORT e AMREX_LANG_FORT
   #
   target_compile_definitions ( amrex PUBLIC
      $<$<COMPILE_LANGUAGE:Fortran>:BL_LANG_FORT> )
   target_compile_definitions ( amrex PUBLIC
      $<$<COMPILE_LANGUAGE:Fortran>:AMREX_LANG_FORT> )

   #
   # GNU-specific defines
   # 
   if ( ${CMAKE_C_COMPILER_ID} STREQUAL "GNU" ) 
   
      if ( CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8" )
	 message ( WARNING " Your default GCC is version ${CMAKE_CXX_COMPILER_VERSION}.This might break during build. GCC>=4.8 is recommended.")
      endif ()
      
      string ( REPLACE "." ";" VERSION_LIST ${CMAKE_CXX_COMPILER_VERSION})
      list ( GET VERSION_LIST 0 GCC_VERSION_MAJOR )
      list ( GET VERSION_LIST 1 GCC_VERSION_MINOR )

      target_compile_definitions ( amrex PUBLIC
	 BL_GCC_VERSION=${CMAKE_CXX_COMPILER_VERSION}
	 BL_GCC_MAJOR_VERSION=${GCC_VERSION_MAJOR}
	 BL_GCC_MINOR_VERSION=${GCC_VERSION_MINOR}
	 )

   endif ()

   # 
   # Fortran/C mangling scheme
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

   add_amrex_define( AMREX_FORT_USE_${FORTLINK} )

   # SENSEI Insitu
   add_amrex_define( AMREX_USE_SENSEI_INSITU IF ENABLE_SENSEI_INSITU )
   
endfunction ()

# 
#
# FUNCTION: add_amrex_define
#
# Add definition to target "amrex" compile definitions.
#
# Arguments:
#
#    new_define    = variable containing the definition to add.
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
# Author: Michele Rosso
# Date  : June 26, 2018
# 
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



