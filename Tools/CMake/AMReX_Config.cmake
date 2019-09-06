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
   include( AMReXGenexHelpers )
   
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
   

   # Set C++ standard and disable compiler-specific extensions, like "-std=gnu++14" for GNU
   # This will also enforce the same standard with the CUDA compiler
   # Moreover, it will also enforce such standard on all the consuming targets
   target_compile_features(amrex PUBLIC cxx_std_14)
   set_target_properties(amrex PROPERTIES CXX_EXTENSIONS OFF) # This disable C++ standard extension
  
   #
   # Setup OpenMP 
   #
   if (ENABLE_OMP)
      find_package(OpenMP REQUIRED CXX Fortran)
      
      target_link_libraries(amrex
         PUBLIC
         OpenMP::OpenMP_CXX
         OpenMP::OpenMP_Fortran)
      
      set_target_properties(
         OpenMP::OpenMP_CXX OpenMP::OpenMP_Fortran
         PROPERTIES
         IMPORTED_GLOBAL True )
      
      # We have to manually pass OpenMP flags to host compiler if CUDA is enabled
      # Since OpenMP imported targets are generated only for the Compiler ID in use, i.e.
      # they do not provide flags for all possible compiler ids, we assume the same compiler used
      # for building amrex will be used to build the application code
      if (ENABLE_CUDA)
         get_target_property(_cxx_omp_flags OpenMP::OpenMP_CXX INTERFACE_COMPILE_OPTIONS)
         
         evaluate_genex(_cxx_omp_flags _omp_flags
            LANG   CXX
            COMP   ${_comp}
            STRING )
         
         target_compile_options(amrex PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${_omp_flags}>)
      endif ()
      
   else ()
      target_compile_options( amrex
         PUBLIC
         $<$<CXX_COMPILER_ID:Cray>:-h;noomp> )     
   endif ()

   
   if (ENABLE_CUDA)     
      # After we load the setups for ALL the supported compilers
      # we can load the setup for NVCC if required
      # This is necessary because we need to know the C++ flags
      # to pass to the Xcompiler option.
      target_compile_definitions( amrex
         PUBLIC
         AMREX_NVCC_VERSION=${CMAKE_CUDA_COMPILER_VERSION}
         AMREX_NVCC_MAJOR_VERSION=${NVCC_VERSION_MAJOR}
         AMREX_NVCC_MINOR_VERSION=${NVCC_VERSION_MINOR} )
      
      set_target_properties( amrex
         PROPERTIES
         CUDA_SEPARABLE_COMPILATION ON      # This adds -dc
         CUDA_RESOLVE_DEVICE_SYMBOLS OFF
         )
      
      #
      # Retrieve compile flags for the current configuration
      # I haven't find a way to set host compiler flags for all the
      # possible configurations.
      #
      get_target_property( _amrex_flags_1 amrex COMPILE_OPTIONS)

      if (NOT CMAKE_CXX_FLAGS)
         get_target_property( _amrex_flags_2 Flags_CXX INTERFACE_COMPILE_OPTIONS)
      endif ()
      
      set(_amrex_flags)
      if (_amrex_flags_1)
         list(APPEND _amrex_flags ${_amrex_flags_1})
      endif ()
      if (_amrex_flags_2)
         list(APPEND _amrex_flags ${_amrex_flags_2})
      endif ()
     
      evaluate_genex(_amrex_flags _amrex_cxx_flags
         LANG   CXX
         COMP   ${CMAKE_CXX_COMPILER_ID}
         CONFIG ${CMAKE_BUILD_TYPE}
         STRING )

      if (_amrex_cxx_flags)
         target_compile_options(amrex PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${_amrex_cxx_flags}>)
      endif ()
      
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
   # Setup third-party profilers
   # 
   set_amrex_profilers()

   #
   # If CUDA is enabled, add manually libcuda because CMake does not find it
   # Do the same for nvToolsExt if tiny profiler is on
   #
   if (ENABLE_CUDA)
      target_link_libraries(amrex PUBLIC cuda)
      if (ENABLE_TINY_PROFILE)
          target_link_libraries(amrex PUBLIC nvToolsExt)
      endif ()
   endif ()

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


  # Include AMReX cmake functions 
  include(AMReXGenexHelpers)
  include(AMReXTargetHelpers)

  get_target_properties_flattened(amrex  _includes _defines _flags _link_line)

  set(_lang CXX Fortran) 
  set(_prop _includes _defines _flags _link_line )


  # Loop over all combinations of language and property and extract 
  # what you need 
  foreach( _p IN LISTS _prop )
     foreach( _l IN LISTS _lang )

        string(TOLOWER ${_l} _ll) # Lower case language name

        # _${_ll}${_p} is a variable named as _lang_property,
        # both lower case. 
        evaluate_genex(${_p} _${_ll}${_p}
           LANG ${_l}
           COMP ${CMAKE_${_l}_COMPILER_ID}
           CONFIG ${CMAKE_BUILD_TYPE}
           INTERFACE BUILD)

        if (_${_ll}${_p})

           list(REMOVE_DUPLICATES _${_ll}${_p})
           
           if ("${_p}" STREQUAL "_defines")
              string(REPLACE ";" " -D" _${_ll}${_p} "-D${_${_ll}${_p}}")
           elseif ("${_p}" STREQUAL "_includes")
              string(REPLACE ";" " -I" _${_ll}${_p} "-I${_${_ll}${_p}}")
           else()
              string(REPLACE ";" " " _${_ll}${_p} "${_${_ll}${_p}}")
           endif ()              

        endif ()
        
     endforeach()
  endforeach ()

   
   string ( TOUPPER "${CMAKE_BUILD_TYPE}"  AMREX_BUILD_TYPE )
   set(_cxx_flags "${CMAKE_CXX_FLAGS_${AMREX_BUILD_TYPE}} ${CMAKE_CXX_FLAGS} ${_cxx_flags}")
   set(_fortran_flags "${CMAKE_Fortran_FLAGS_${AMREX_BUILD_TYPE}} ${CMAKE_Fortran_FLAGS} ${_fortran_flags}")   

   
   #
   # Config summary
   #
   message( STATUS "AMReX configuration summary: ")
   message( STATUS "   Build type               = ${CMAKE_BUILD_TYPE}")
   message( STATUS "   Install directory        = ${CMAKE_INSTALL_PREFIX}")
   message( STATUS "   C++ compiler             = ${CMAKE_CXX_COMPILER}")
   message( STATUS "   Fortran compiler         = ${CMAKE_Fortran_COMPILER}")
   if (ENABLE_CUDA)
      message( STATUS "   CUDA compiler            = ${CMAKE_CUDA_COMPILER}")
   endif ()
   
   message( STATUS "   C++ defines              = ${_cxx_defines}")
   message( STATUS "   Fortran defines          = ${_fortran_defines}")
   
   message( STATUS "   C++ flags                = ${_cxx_flags}")
   message( STATUS "   Fortran flags            = ${_fortran_flags}")
   if (ENABLE_CUDA)
      message( STATUS "   CUDA flags               = ${CMAKE_CUDA_FLAGS_${AMREX_BUILD_TYPE}} ${CMAKE_CUDA_FLAGS}"
         "${AMREX_CUDA_FLAGS}")
   endif ()
   message( STATUS "   C++ include paths        = ${_cxx_includes}")  
   message( STATUS "   Fortran include paths    = ${_fortran_includes}")
   message( STATUS "   Link line                = ${_cxx_link_line}") 

endfunction ()

