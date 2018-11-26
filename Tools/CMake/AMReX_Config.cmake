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
   # Setup MPI (CMake 3.11 allows to import MPI as an imported TARGET)
   #  
   if (ENABLE_MPI)
      find_package (MPI REQUIRED)

      # Includes
      target_include_directories ( amrex PUBLIC
	 ${MPI_Fortran_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH}
	 ${MPI_CXX_INCLUDE_PATH} )

      # Genex $<COMPILE_LANGUAGE:...> is broken for export. A fix
      # is present in CMake 3.12
      # target_include_directories ( amrex PUBLIC
      # 	 $<$<COMPILE_LANGUAGE:Fortran>:${MPI_Fortran_INCLUDE_PATH}>
      # 	 $<$<COMPILE_LANGUAGE:C>:${MPI_C_INCLUDE_PATH}>
      # 	 $<$<COMPILE_LANGUAGE:CXX>:${MPI_CXX_INCLUDE_PATH}>
      # 	 )

      
      # Additional compiler flags
      strip (MPI_C_COMPILE_FLAGS)
      strip (MPI_CXX_COMPILE_FLAGS)
      strip (MPI_Fortran_COMPILE_FLAGS)
      target_compile_options ( amrex PUBLIC
	 $<$<COMPILE_LANGUAGE:Fortran>:${MPI_Fortran_COMPILE_FLAGS}>
	 $<$<COMPILE_LANGUAGE:C>:${MPI_C_COMPILE_FLAGS}>
	 $<$<COMPILE_LANGUAGE:CXX>:${MPI_CXX_COMPILE_FLAGS}>
	 )
      
      # Libraries to link to + linker flags
      # CMake doc suggests not to use absolute paths for portability reasons
      # (and we ignore it for now)
      strip ( MPI_CXX_LINK_FLAGS )
      target_link_libraries ( amrex PUBLIC ${MPI_CXX_LINK_FLAGS}
	 ${MPI_Fortran_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
      
   endif ()

   #
   # Setup OpenMP
   #
   if (ENABLE_OMP)
      find_package (OpenMP REQUIRED)

      # Additional compiler flags
      strip (OpenMP_C_FLAGS)
      strip (OpenMP_CXX_FLAGS)
      strip (OpenMP_Fortran_FLAGS)
      target_compile_options ( amrex PUBLIC
	 $<$<COMPILE_LANGUAGE:Fortran>:${OpenMP_Fortran_FLAGS}>
	 $<$<COMPILE_LANGUAGE:C>:${OpenMP_C_FLAGS}>
	 $<$<COMPILE_LANGUAGE:CXX>:${OpenMP_CXX_FLAGS}>
	 )

      # The openMP flags are required also during linking phase
      # The NON-clean way of achiving this is by adding a flag to the linker
      # as follows. Notice that in more recent version of CMake, openMP becomes
      # an imported target, thus provinding a clean way to do this
      strip (OpenMP_CXX_FLAGS)
      target_link_libraries ( amrex PUBLIC ${OpenMP_CXX_FLAGS} )
   else ()
      # Cray compiler has OMP turned on by default
      target_compile_options ( amrex PUBLIC $<$<CXX_COMPILER_ID:Cray>:-h;noomp> $<$<C_COMPILER_ID:Cray>:-h;noomp> )
   endif()
      
   #
   # Add third party libraries
   #
   if (ENABLE_3D_NODAL_MLMG)
      # Add Algoim dependency
      target_include_directories( amrex PRIVATE ${ALGOIM_INSTALL_DIR}/src )
      #  Blitz dependency
      target_include_directories( amrex PRIVATE ${BLITZ_INSTALL_DIR}/include )
      target_link_libraries( amrex PUBLIC ${BLITZ_INSTALL_DIR}/lib/libblitz.a)
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
   message( STATUS "   C++ flags                = ${AMREX_CXX_FLAGS}")
   message( STATUS "   Fortran flags            = ${AMREX_Fortran_FLAGS}")
   message( STATUS "   C++ include paths        = ${AMREX_CXX_INCLUDE_PATHS}")  
   message( STATUS "   Fortran include paths    = ${AMREX_Fortran_INCLUDE_PATHS}")
   message( STATUS "   Link line                = ${AMREX_LINK_LINE}") 

endfunction ()
