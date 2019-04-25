#
#
# FUNCTION: configure_amrex
#
# Set all the properties (except sources and installation-related)
# required to build target "amrex".
# Target "amrex" must exist before this function is called.
#
# Author: Michele Rosso
# Date  : June 26, 2018
#
#
function (configure_amrex)

   # 
   # Check if target "amrex" has been defined before
   # calling this macro
   #
   if (NOT TARGET amrex)
      message (FATAL_ERROR "Target 'amrex' must be defined before calling function 'configure_amrex'" )     
   endif ()

   #
   # Check that needed options have already been defined
   # 
   if ( ( NOT ( DEFINED ENABLE_MPI ) ) OR ( NOT (DEFINED ENABLE_OMP) ) 
	 OR ( NOT (DEFINED ENABLE_PIC) ) OR (NOT (DEFINED ENABLE_FPE)))
      message ( AUTHOR_WARNING "Required options are not defined" )
   endif ()

   #
   # Include the required modules
   # 
   include( AMReX_ThirdPartyProfilers )
   include( AMReX_Defines )
   include( AMReX_Utils )
   
   # 
   # Set properties for target "amrex"
   # 
   set_amrex_defines()
   
   #
   # Setup compilers
   #
  
   # Exit if Cray compiler is in use -- Support for Cray is currently broken   
   if ( ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray") OR
         ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Cray") )
      message(FATAL_ERROR "Support for Cray compiler is currently broken")
   endif()

   # Load flags preset
   set_compiler_flags_preset(amrex)

   if (ENABLE_CUDA)
      # After we load the setups for ALL the supported compilers
      # we can load the setup for NVCC if required
      # This is necessary because we need to know the C++ flags
      # to pass to the Xcompiler option.
      # include(${CMAKE_MODULE_PATH}/comps/AMReX_NVIDIA.cmake)
      target_compile_definitions( amrex PUBLIC
         AMREX_NVCC_VERSION=${CMAKE_CUDA_COMPILER_VERSION}
         AMREX_NVCC_MAJOR_VERSION=${NVCC_VERSION_MAJOR}
         AMREX_NVCC_MINOR_VERSION=${NVCC_VERSION_MINOR} )
      
      set_target_properties( amrex
         PROPERTIES
         CUDA_SEPARABLE_COMPILATION ON      # This adds -dc
         CUDA_RESOLVE_DEVICE_SYMBOLS OFF
         )
      
      if (NOT ENABLE_3D_NODAL_MLMG)
         set_target_properties( amrex
            PROPERTIES
            CUDA_STANDARD 11     # Adds -std=<standard>
            CUDA_STANDARD_REQUIRED ON
            )
      else ()
         set_target_properties( amrex
            PROPERTIES
            CUDA_STANDARD 14     # Adds -std=<standard>
            CUDA_STANDARD_REQUIRED ON
            )
      endif ()
      
      #
      # Retrieve compile flags for the current configuration
      # I haven't find a way to set host compiler flags for all the
      # possible configurations.
      #
      get_target_property( _amrex_flags amrex COMPILE_OPTIONS)

      if (NOT CMAKE_CXX_FLAGS)
         get_target_property( _amrex_flags_2 Flags_CXX INTERFACE_COMPILE_OPTIONS)
         list(APPEND _amrex_flags ${_amrex_flags_2})
      endif ()

      evaluate_genex(_amrex_flags _amrex_cxx_flags
         LANG   CXX
         COMP   ${CMAKE_CXX_COMPILER_ID}
         CONFIG ${CMAKE_BUILD_TYPE}
         STRING )

      target_compile_options(amrex PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${_amrex_cxx_flags}>)

   endif ()
   

   #
   # GNU-specific defines
   # 
   if ( ${CMAKE_C_COMPILER_ID} STREQUAL "GNU" ) 
      
      if ( CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8" )
         message( WARNING
            " Your default GCC is version ${CMAKE_CXX_COMPILER_VERSION}.This might break during build. GCC>=4.8 is recommended.")
      endif ()
      
      string( REPLACE "." ";" VERSION_LIST ${CMAKE_CXX_COMPILER_VERSION})
      list( GET VERSION_LIST 0 GCC_VERSION_MAJOR )
      list( GET VERSION_LIST 1 GCC_VERSION_MINOR )

      target_compile_definitions( amrex PUBLIC $<BUILD_INTERFACE:
         BL_GCC_VERSION=${CMAKE_CXX_COMPILER_VERSION}
         BL_GCC_MAJOR_VERSION=${GCC_VERSION_MAJOR}
         BL_GCC_MINOR_VERSION=${GCC_VERSION_MINOR}
         >
         )
   endif ()


   if ( ENABLE_PIC OR BUILD_SHARED_LIBS )
      set_target_properties ( amrex PROPERTIES POSITION_INDEPENDENT_CODE True )
   endif ()

   if ( BUILD_SHARED_LIBS OR ENABLE_CUDA )
      if(APPLE)
         target_link_options(amrex PUBLIC -Wl,-undefined,warning)
      else()
         target_link_options(amrex PUBLIC -Wl,--warn-unresolved-symbols)
      endif()
   endif() 
   
   #
   # Setup OpenMP 
   #
   if (ENABLE_OMP)
      find_package(OpenMP REQUIRED)
      target_link_libraries(amrex PUBLIC OpenMP::OpenMP_CXX)
   else ()
      target_compile_options( amrex PUBLIC
         $<$<CXX_COMPILER_ID:Cray>:-h;noomp> )     
   endif ()
   
   #
   # Add third party libraries
   #
   if (ENABLE_3D_NODAL_MLMG)
      include(AMReX_InstallExternalLibs) 
   endif()
   
   #
   # Setup third-party profilers
   # 
   set_amrex_profilers()

   #
   # If CUDA is enabled, add manually libcuda because CMake does not find it
   #
   if (ENABLE_CUDA)
      target_link_libraries(amrex PUBLIC cuda)   
   endif ()

   #
   # Print out summary
   #  
   print_amrex_configuration_summary ()

endfunction ()

# 
#
# Prints out configuration details
# 
# 
function (print_amrex_configuration_summary)

   # 
   # Check if target "amrex" has been defined before
   # calling this macro
   #
   if (NOT TARGET amrex)
      message (FATAL_ERROR "Target 'amrex' must be defined before calling function 'set_amrex_defines'" )
   endif ()
   
   #
   # Retrieve defines
   #
   get_target_property ( AMREX_DEFINES amrex COMPILE_DEFINITIONS )
   evaluate_genex(AMREX_DEFINES AMREX_CXX_DEFINES     LANG CXX     COMP ${CMAKE_CXX_COMPILER_ID} )
   evaluate_genex(AMREX_DEFINES AMREX_Fortran_DEFINES LANG Fortran COMP ${CMAKE_Fortran_COMPILER_ID} ) 
   string (REPLACE ";" " -D" AMREX_CXX_DEFINES "-D${AMREX_CXX_DEFINES}" )
   string (REPLACE ";" " -D" AMREX_Fortran_DEFINES "-D${AMREX_Fortran_DEFINES}" )
   
   #
   # Retrieve compile flags
   #
   get_target_property( AMREX_FLAGS   amrex COMPILE_OPTIONS)

   evaluate_genex(AMREX_FLAGS AMREX_CXX_FLAGS
      LANG   CXX
      COMP   ${CMAKE_CXX_COMPILER_ID}
      CONFIG ${CMAKE_BUILD_TYPE}
      STRING )

   evaluate_genex(AMREX_FLAGS AMREX_Fortran_FLAGS
      LANG   Fortran
      COMP   ${CMAKE_Fortran_COMPILER_ID}
      CONFIG ${CMAKE_BUILD_TYPE}
      STRING )

   if (ENABLE_CUDA)
      evaluate_genex(AMREX_FLAGS AMREX_CUDA_FLAGS
         LANG   CUDA
         COMP   ${CMAKE_CUDA_COMPILER_ID}
         CONFIG ${CMAKE_BUILD_TYPE}
         STRING )
   endif ()
   
   # Add base flags
   string ( TOUPPER "${CMAKE_BUILD_TYPE}"  AMREX_BUILD_TYPE )
   set (AMREX_CXX_FLAGS  "${CMAKE_CXX_FLAGS_${AMREX_BUILD_TYPE}} ${CMAKE_CXX_FLAGS} ${AMREX_CXX_FLAGS}")
   set (AMREX_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_${AMREX_BUILD_TYPE}} ${CMAKE_Fortran_FLAGS} ${AMREX_Fortran_FLAGS}")
   string (STRIP "${AMREX_CXX_FLAGS}" AMREX_CXX_FLAGS)
   string (STRIP "${AMREX_Fortran_FLAGS}" AMREX_Fortran_FLAGS)

   #
   # Include paths
   #
   get_target_property( AMREX_INCLUDE_PATHS amrex INTERFACE_INCLUDE_DIRECTORIES )
   evaluate_genex(AMREX_INCLUDE_PATHS AMREX_CXX_INCLUDE_PATHS     LANG CXX )
   evaluate_genex(AMREX_INCLUDE_PATHS AMREX_Fortran_INCLUDE_PATHS LANG Fortran )

   #
   # Link libraries
   # 
   get_target_property ( TMP amrex LINK_LIBRARIES )
   evaluate_genex(TMP AMREX_LINK_LINE )
   if (NOT AMREX_LINK_LINE) # LINK_LIBRARIES property can return "NOT_FOUND"
      set (AMREX_LINK_LINE "")
   endif ()   
   string ( REPLACE ";" " " AMREX_LINK_LINE "${AMREX_LINK_LINE}" )

   
   #
   # Config summary
   #
   message( STATUS "AMReX configuration summary: ")
   message( STATUS "   Build type               = ${CMAKE_BUILD_TYPE}")
   message( STATUS "   Install directory        = ${CMAKE_INSTALL_PREFIX}")
   message( STATUS "   C++ defines              = ${AMREX_CXX_DEFINES}")
   message( STATUS "   Fortran defines          = ${AMREX_Fortran_DEFINES}")
   message( STATUS "   C++ compiler             = ${CMAKE_CXX_COMPILER}")
   message( STATUS "   Fortran compiler         = ${CMAKE_Fortran_COMPILER}")
   if (ENABLE_CUDA)
      message( STATUS "   CUDA compiler            = ${CMAKE_CUDA_COMPILER}")
   endif ()
   message( STATUS "   C++ flags                = ${AMREX_CXX_FLAGS}")
   message( STATUS "   Fortran flags            = ${AMREX_Fortran_FLAGS}")
   if (ENABLE_CUDA)
      message( STATUS "   CUDA flags               = ${CMAKE_CUDA_FLAGS_${AMREX_BUILD_TYPE}} ${CMAKE_CUDA_FLAGS}"
         "${AMREX_CUDA_FLAGS}")
   endif ()
   message( STATUS "   C++ include paths        = ${AMREX_CXX_INCLUDE_PATHS}")  
   message( STATUS "   Fortran include paths    = ${AMREX_Fortran_INCLUDE_PATHS}")
   message( STATUS "   Link line                = ${AMREX_LINK_LINE}") 

endfunction ()


#
# For now let's keep this here
# 
function ( set_compiler_flags_preset _target )

   #
   # WARNING: since cmake 3.14 the genex Fortran_COMPILER_ID is available!!!
   #

   # 
   # Check if target "_target" has been defined before
   # calling this macro
   #
   if ( NOT TARGET ${_target} )
      message (FATAL_ERROR "Target '${_target}' does not exist" )
   endif ()

   #
   # Check wether the compiler ID has been defined
   # 
   if (  NOT (DEFINED CMAKE_Fortran_COMPILER_ID) OR
	 NOT (DEFINED CMAKE_C_COMPILER_ID) OR 
	 NOT (DEFINED CMAKE_CXX_COMPILER_ID) )
      message ( FATAL_ERROR "Compiler ID is UNDEFINED" )
   endif ()

   #
   # Helper variables
   # 
   set(_cxx             "$<COMPILE_LANGUAGE:CXX>")
   set(_fortran         "$<COMPILE_LANGUAGE:Fortran>")

   set(_debug           "$<CONFIG:Debug>")
   set(_release         "$<CONFIG:Release>")

   set(_gnu             "$<CXX_COMPILER_ID:GNU>")
   set(_cxx_gnu         "$<AND:${_cxx},${_gnu}>")

   set(_intel           "$<CXX_COMPILER_ID:Intel>")
   set(_cxx_intel       "$<AND:${_cxx},${_intel}>")

   set(_pgi             "$<CXX_COMPILER_ID:PGI>")
   set(_cxx_pgi         "$<AND:${_cxx},${_pgi}>")

   set(_cray            "$<CXX_COMPILER_ID:Cray>")
   set(_cxx_cray        "$<AND:${_cxx},${_cray}>")

   set(_clang           "$<CXX_COMPILER_ID:Clang>")
   set(_cxx_clang       "$<AND:${_cxx},${_clang}>")
   
   set(_apple           "$<CXX_COMPILER_ID:AppleClang>")
   set(_cxx_apple       "$<AND:${_cxx},${_apple}>")


   if (ENABLE_3D_NODAL_MLMG)
      set(_cxx_std c++14)
   else ()
      set(_cxx_std c++11)
   endif ()
   
   target_compile_options(  ${_target}
         PUBLIC
         $<${_cxx_gnu}:-std=${_cxx_std}>
         $<${_cxx_intel}:-std=${_cxx_std}>
         $<${_cxx_cray}:-h std=${_cxx_std} -h list=a>
         $<${_cxx_pgi}:-std=${_cxx_std}>
         $<${_cxx_clang}:-std=${_cxx_std}>
         $<${_cxx_apple}:-std=${_cxx_std}>
         )
  
endfunction () 
