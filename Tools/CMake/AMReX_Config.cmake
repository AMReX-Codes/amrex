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
	 OR ( NOT (DEFINED ENABLE_PIC ) ) )
      message ( AUTHOR_WARNING "Required options are not defined" )
   endif ()

   #
   # Include the required modules
   # 
   include( AMReX_ThirdPartyProfilers )
   include( AMReX_Defines )
   include( AMReX_Compilers )
   include( AMReX_Utils )
   
   # 
   # Set properties for target "amrex"
   # 
   set_amrex_defines ()
   set_amrex_compiler_flags ()

   if ( ENABLE_PIC OR BUILD_SHARED_LIBS )
      set_target_properties ( amrex PROPERTIES POSITION_INDEPENDENT_CODE True )
   endif ()
      
   # 
   # Location of Fortran modules
   # 
   set ( AMREX_Fortran_MODULE_DIR ${PROJECT_BINARY_DIR}/mod_files )

   set_target_properties ( amrex
      PROPERTIES
      Fortran_MODULE_DIRECTORY ${AMREX_Fortran_MODULE_DIR} )
   
   target_include_directories ( amrex
      PUBLIC "$<BUILD_INTERFACE:${AMREX_Fortran_MODULE_DIR}>" )
   
   #
   # Setup MPI (CMake >= 3.11 allows to import MPI as an imported TARGET)
   # Check https://cliutils.gitlab.io/modern-cmake/chapters/packages/MPI.html
   #  
   if (ENABLE_MPI)
      find_package(MPI REQUIRED)
      target_link_libraries(amrex PUBLIC MPI::MPI_CXX MPI::MPI_Fortran)   
   endif ()
   
   #
   # Setup OpenMP (CMake >= 3.9 allows to import MPI as an imported TARGET)
   #
   if (ENABLE_OMP)
      find_package(OpenMP REQUIRED)

      #
      # This is to allow building on macOS -- This is not needed with CMake 3.12+
      # Check https://cliutils.gitlab.io/modern-cmake/chapters/packages/OpenMP.html
      #
      if(${CMAKE_MINIMUM_REQUIRED_VERSION} VERSION_GREATER_EQUAL "3.12.0") 
         message(AUTHOR_WARNING "Remove CMAke legacy support for OpenMP")
      else()
         if(NOT TARGET OpenMP::OpenMP_CXX)
            find_package(Threads REQUIRED)
            add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
            set_property(TARGET OpenMP::OpenMP_CXX
               PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
            # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
            set_property(TARGET OpenMP::OpenMP_CXX
               PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)
         endif()       
      endif()

      # I don't think there's a specific Fortran library to link for pthreads
      target_link_libraries(amrex PUBLIC OpenMP::OpenMP_CXX)

   else ()
      # Cray compiler has OMP turned on by default
      target_compile_options( amrex PUBLIC
         $<$<CXX_COMPILER_ID:Cray>:-h;noomp> $<$<C_COMPILER_ID:Cray>:-h;noomp> )     
   endif ()

   #
   # Set compile features -- We just need to set C++ standard
   #
   set(SUPPORTED_COMPILERS GNU Intel PGI)
   
   if ("${CMAKE_CXX_COMPILER_ID}" IN_LIST SUPPORTED_COMPILERS)
      if (ENABLE_3D_NODAL_MLMG)
         target_compile_features(amrex PUBLIC cxx_std_14)
      else ()
         target_compile_features(amrex PUBLIC cxx_std_11)
      endif ()
   else ()
      # This branch is for compiler not yet supported by CMake compiler features
      # Remove this branch as more compilers are supported
      if (NOT ENABLE_3D_NODAL_MLMG)
         target_compile_options ( amrex
            PUBLIC
            $<$<CXX_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:CXX>:-h std=c++11 -h list=a>>
            $<$<CXX_COMPILER_ID:Clang>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>>
            $<$<CXX_COMPILER_ID:AppleClang>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>> )
      else ()
         target_compile_options ( amrex
            PUBLIC
            $<$<CXX_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:CXX>:-h std=c++14 -h list=a>>
            $<$<CXX_COMPILER_ID:Clang>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++14>>
            $<$<CXX_COMPILER_ID:AppleClang>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++14>> )
      endif ()
   endif ()
   
   #
   # Add third party libraries
   #
   if (ENABLE_3D_NODAL_MLMG)
      include(AMReX_InstallExternalLibs)
      target_link_libraries(amrex PRIVATE algoim PUBLIC blitz)
   endif()
   
   #
   # Setup third-party profilers
   # 
   set_amrex_profilers()

   #
   # Setup SENSEI
   #
   if (ENABLE_SENSEI_INSITU)
      find_package(SENSEI REQUIRED)
   endif()

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
   replace_genex ( AMREX_DEFINES AMREX_CXX_DEFINES     LANGUAGE CXX )
   replace_genex ( AMREX_DEFINES AMREX_Fortran_DEFINES LANGUAGE Fortran )


   # extract_by_language ( AMREX_DEFINES AMREX_CXX_DEFINES AMREX_Fortran_DEFINES )
   string (REPLACE " " " -D" AMREX_CXX_DEFINES "-D${AMREX_CXX_DEFINES}" )
   string (REPLACE ";" " -D" AMREX_CXX_DEFINES "${AMREX_CXX_DEFINES}" )
   string (REPLACE " " " -D" AMREX_Fortran_DEFINES "-D${AMREX_Fortran_DEFINES}" )
   string (REPLACE ";" " -D" AMREX_Fortran_DEFINES "${AMREX_Fortran_DEFINES}" )
   
   #
   # Retrieve compile flags
   #
   get_target_property ( AMREX_FLAGS   amrex COMPILE_OPTIONS)
   replace_genex ( AMREX_FLAGS AMREX_CXX_FLAGS     LANGUAGE CXX )
   replace_genex ( AMREX_FLAGS AMREX_Fortran_FLAGS LANGUAGE Fortran )
   string ( REPLACE ";" " " AMREX_CXX_FLAGS "${AMREX_CXX_FLAGS}" )
   string ( REPLACE ";" " " AMREX_Fortran_FLAGS "${AMREX_Fortran_FLAGS}" )

   # Add base flags
   string ( TOUPPER "${CMAKE_BUILD_TYPE}"  AMREX_BUILD_TYPE )
   set (AMREX_CXX_FLAGS  "${CMAKE_CXX_FLAGS_${AMREX_BUILD_TYPE}} ${CMAKE_CXX_FLAGS} ${AMREX_CXX_FLAGS}")
   set (AMREX_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_${AMREX_BUILD_TYPE}} ${CMAKE_Fortran_FLAGS} ${AMREX_Fortran_FLAGS}")
   string (STRIP "${AMREX_CXX_FLAGS}" AMREX_CXX_FLAGS)
   string (STRIP "${AMREX_Fortran_FLAGS}" AMREX_Fortran_FLAGS)

   #
   # Include paths
   #
   get_target_property ( AMREX_INCLUDE_PATHS amrex INTERFACE_INCLUDE_DIRECTORIES )
   replace_genex ( AMREX_INCLUDE_PATHS AMREX_CXX_INCLUDE_PATHS LANGUAGE CXX )
   replace_genex ( AMREX_INCLUDE_PATHS AMREX_Fortran_INCLUDE_PATHS LANGUAGE Fortran )


   #
   # Link libraries
   # 
   get_target_property ( TMP amrex LINK_LIBRARIES )
   replace_genex (TMP AMREX_LINK_LINE)
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
      message( STATUS "   CUDA flags               = ${CMAKE_CUDA_FLAGS_${AMREX_BUILD_TYPE}} ${CMAKE_CUDA_FLAGS}")
   endif ()
   message( STATUS "   C++ include paths        = ${AMREX_CXX_INCLUDE_PATHS}")  
   message( STATUS "   Fortran include paths    = ${AMREX_Fortran_INCLUDE_PATHS}")
   message( STATUS "   Link line                = ${AMREX_LINK_LINE}") 

endfunction ()
