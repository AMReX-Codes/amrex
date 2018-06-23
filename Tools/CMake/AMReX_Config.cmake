###############################################

# Here we configure the build                 #

###############################################
if ( DEFINED __AMREX_CONFIG__ )
   return ()
endif ()
set ( __AMREX_CONFIG__ )

#
# Check if AMReX_Options.cmake  and AMREX_CMakeVariables.cmake
# have been already processed
#
if ( NOT (DEFINED __AMREX_OPTIONS__ ) )
   message ( FATAL_ERROR "AMReX_Options.cmake must be\
included before AMReX_Configure.cmake" )
endif ()

if ( NOT (DEFINED __AMREX_CMAKEVARIABLES__ ) )
   message ( FATAL_ERROR "AMReX_CMakeVariables.cmake must be\
included before AMReX_Configure.cmake" )
endif ()

#
# Decide whether or not to use PIC 
#
if ( ENABLE_PIC OR BUILD_SHARED_LIBS )
   set (CMAKE_POSITION_INDEPENDENT_CODE TRUE)
endif ()

#
# Setup third-party profilers
#
include ( AMReX_ThirdPartyProfilers )

# Includes
list (APPEND AMREX_EXTRA_Fortran_INCLUDE_PATH "${TPP_Fortran_INCLUDE_PATH}")
list (APPEND AMREX_EXTRA_C_INCLUDE_PATH "${TPP_C_INCLUDE_PATH}")
list (APPEND AMREX_EXTRA_CXX_INCLUDE_PATH "${TPP_CXX_INCLUDE_PATH}")

# Compile flags
append ( TPP_FFLAGS AMREX_EXTRA_Fortran_FLAGS ) 
append ( TPP_CFLAGS AMREX_EXTRA_C_FLAGS )
append ( TPP_CXXFLAGS AMREX_EXTRA_CXX_FLAGS )

# Link Line
append_to_link_line ( TPP_Fortran_LINK_LINE AMREX_EXTRA_Fortran_LINK_LINE )
append_to_link_line ( TPP_C_LINK_LINE AMREX_EXTRA_C_LINK_LINE )
append_to_link_line ( TPP_CXX_LINK_LINE AMREX_EXTRA_CXX_LINK_LINE )

#
# Set preprocessor flags
#
include ( AMReX_Defines )

#
# Setup compiler flags 
#
include ( AMReX_Compilers )

# If CMAKE_<LANG>_FLAGS is not defined by user, load defaults
if ( NOT CMAKE_CXX_FLAGS )
   set ( CMAKE_CXX_FLAGS "${AMREX_CXXFLAGS_${AMREX_BUILD_TYPE}}" )
endif ()

if ( NOT CMAKE_Fortran_FLAGS )
   set ( CMAKE_Fortran_FLAGS "${AMREX_FFLAGS_${AMREX_BUILD_TYPE}}" )
endif ()

# Use CMAKE_<LANG>_FLAGS_<BUILD_TYPE> to store required flags
# NOTE: this is not the standard CMake way of doing things
set ( CMAKE_CXX_FLAGS_${AMREX_BUILD_TYPE}  "${AMREX_CXXFLAGS_REQUIRED}" )
set ( CMAKE_Fortran_FLAGS_${AMREX_BUILD_TYPE}  "${AMREX_FFLAGS_REQUIRED}" )

if ( ENABLE_FPE )
   append ( AMREX_CXXFLAGS_FPE CMAKE_CXX_FLAGS_${AMREX_BUILD_TYPE} )
   append ( AMREX_FFLAGS_FPE   CMAKE_Fortran_FLAGS_${AMREX_BUILD_TYPE} )
endif ()

#  Add definition related to specific compiler ( only GNU for now )
#  AMREX_COMPILER_DEFINES is defined in AMReX_Compilers.cmake
if ( AMREX_COMPILER_DEFINES )
   foreach ( item ${AMREX_COMPILER_DEFINES} )
      add_define (${item})
   endforeach ()
endif ()


# Add extra flags
append ( AMREX_EXTRA_Fortran_FLAGS CMAKE_Fortran_FLAGS )
append ( AMREX_EXTRA_CXX_FLAGS CMAKE_CXX_FLAGS )

# Accumulate all the flags into AMREX_<LANG>_FLAGS: these variables will be exported
set ( AMREX_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${AMREX_BUILD_TYPE}}" )
set ( AMREX_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${AMREX_BUILD_TYPE}}" )

#
# Config summary
#
message( STATUS "AMReX configuration summary: ")
message( STATUS "   Build type               = ${CMAKE_BUILD_TYPE}")
message( STATUS "   Install directory        = ${CMAKE_INSTALL_PREFIX}")
message( STATUS "   Preprocessor flags       = ${AMREX_DEFINES}")
message( STATUS "   C++ compiler             = ${CMAKE_CXX_COMPILER}")
message( STATUS "   Fortran compiler         = ${CMAKE_Fortran_COMPILER}")
message( STATUS "   C++ flags                = ${AMREX_CXX_FLAGS}")
message( STATUS "   Fortran flags            = ${AMREX_Fortran_FLAGS}")
message( STATUS "   C++ include paths        = ${AMREX_EXTRA_CXX_INCLUDE_PATH}") 
message( STATUS "   Fortran include paths    = ${AMREX_EXTRA_Fortran_INCLUDE_PATH}")
message( STATUS "   C++ external libs        = ${AMREX_EXTRA_CXX_LIBRARIES}") 
message( STATUS "   Fortran external libs    = ${AMREX_EXTRA_Fortran_LIBRARIES}")
message( STATUS "   C++ link flags           = ${AMREX_EXTRA_CXX_LINK_FLAGS}") 
message( STATUS "   Fortran link flags       = ${AMREX_EXTRA_Fortran_LINK_FLAGS}")
message( STATUS "   C++ extra link line      = ${AMREX_EXTRA_CXX_LINK_LINE}") 
message( STATUS "   Fortran extra link line  = ${AMREX_EXTRA_Fortran_LINK_LINE}")



function (configure_amrex)

   # 
   # Check if target "amrex" has been defined before
   # calling this macro
   #
   if (NOT TARGET amrex)
      message (FATAL_ERROR "Target 'amrex' must be defined before calling function 'configure_amrex'" )     
   endif ()

   #
   # MPI (CMake 3.11 allows to import MPI as an imported TARGET)
   #  
   if (ENABLE_MPI)
      find_package (MPI REQUIRED)

      # Includes
      target_include_directories ( amrex PUBLIC
	 $<$<COMPILE_LANGUAGE:Fortran>:${MPI_Fortran_INCLUDE_PATH}>
	 $<$<COMPILE_LANGUAGE:C>:${MPI_C_INCLUDE_PATH}>
	 $<$<COMPILE_LANGUAGE:CXX>:${MPI_CXX_INCLUDE_PATH}>
	 )

      # Additional compiler flags
      target_compile_options ( amrex PUBLIC
	 $<$<COMPILE_LANGUAGE:Fortran>:${MPI_Fortran_COMPILE_FLAGS}>
	 $<$<COMPILE_LANGUAGE:C>:${MPI_C_COMPILE_FLAGS}>
	 $<$<COMPILE_LANGUAGE:CXX>:${MPI_CXX_COMPILE_FLAGS}>
	 )
      
      # Libraries to link to + linker flags
      string ( STRIP ${MPI_CXX_LINK_FLAGS} MPI_CXX_LINK_FLAGS )
      target_link_libraries ( amrex PUBLIC "${MPI_CXX_LINK_FLAGS}"
	 ${MPI_Fortran_LIBRARIES} ${MPI_C_LIBRARIES} )
	 
   endif ()


   #
   # Setup OpenMP
   # 
   if (ENABLE_OMP)
      find_package (OpenMP REQUIRED)

      # Additional compiler flags
      target_compile_options ( amrex PUBLIC
	 $<$<COMPILE_LANGUAGE:Fortran>:${OpenMP_Fortran_FLAGS}>
	 $<$<COMPILE_LANGUAGE:C>:${OpenMP_C_FLAGS}>
	 $<$<COMPILE_LANGUAGE:CXX>:${OpenMP_CXX_FLAGS}>
	 )

      # The openMP flags are required also during linking phase
      # The NON-clean way of achiving this is by adding a flag to the linker
      # as follows. Notice that in more recent version of CMake, openMP becomes
      # an imported target, thus provinding a clean way to do this
      string ( STRIP ${OpenMP_CXX_FLAGS} OpenMP_CXX_FLAGS )
      target_link_libraries ( amrex PUBLIC "${OpenMP_CXX_FLAGS}" )

   # else ()
   #    if ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray" ) # Cray has OMP on by default
   # 	 list ( APPEND AMREX_EXTRA_Fortran_FLAGS  "-h noomp")
   #    endif ()
   #    if ( ${CMAKE_C_COMPILER_ID} STREQUAL "Cray" ) # Cray has OMP on by default
   # 	 list ( APPEND AMREX_EXTRA_C_FLAGS  "-h noomp")
   #    endif ()
   #    if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Cray" ) # Cray has OMP on by default
   # 	 list ( APPEND AMREX_EXTRA_CXX_FLAGS  "-h noomp")
   #    endif ()
   endif()



   # 
   # Set compile definitions
   # 
   set_amrex_defines ()

   #
   # Set compiler flags
   #
   set_amrex_compilers ()
   
endfunction ()
