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
   # Set C++ standard and disable compiler-specific extensions, like "-std=gnu++14" for GNU
   # This will also enforce the same standard with the CUDA compiler
   # Moreover, it will also enforce such standard on all the consuming targets
   #
   set_target_properties(amrex PROPERTIES CXX_EXTENSIONS OFF)
   # minimum: C++11 on Linux, C++17 on Windows
   target_compile_features(amrex PUBLIC $<IF:$<STREQUAL:$<PLATFORM_ID>,Windows>,cxx_std_17,cxx_std_11>)

   if (ENABLE_CUDA AND (CMAKE_VERSION VERSION_GREATER_EQUAL 3.17) )
      set_target_properties(amrex PROPERTIES CUDA_EXTENSIONS OFF)
      # minimum: C++11 on Linux, C++17 on Windows
      target_compile_features(amrex PUBLIC $<IF:$<STREQUAL:$<PLATFORM_ID>,Windows>,cuda_std_17,cuda_std_11>)
   endif()

   #
   # Special flags for MSV compiler
   #
   set(_cxx_msvc   "$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:MSVC>>")
   set(_condition  "$<VERSION_LESS:$<CXX_COMPILER_VERSION>,19.26>")

   target_compile_options( amrex PRIVATE $<${_cxx_msvc}:/bigobj> )

   target_compile_options( amrex PUBLIC
      $<${_cxx_msvc}:$<IF:${_condition},/experimental:preprocessor,/Zc:preprocessor>>
      )

   unset(_condition)
   unset(_cxx_msvc)

   #
   # Setup OpenMP
   #
   if (ENABLE_OMP)
      # We have to manually pass OpenMP flags to host compiler if CUDA is enabled
      # Since OpenMP imported targets are generated only for the Compiler ID in use, i.e.
      # they do not provide flags for all possible compiler ids, we assume the same compiler use
      # for building amrex will be used to build the application code
      if (ENABLE_CUDA)
         get_target_property(_omp_flags OpenMP::OpenMP_CXX INTERFACE_COMPILE_OPTIONS)

         eval_genex(_omp_flags CXX ${_comp} INTERFACE BUILD STRING )

         target_compile_options(amrex PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${_omp_flags}>)
      endif ()

   else ()
      target_compile_options( amrex
         PUBLIC
         $<$<CXX_COMPILER_ID:Cray>:-h;noomp> )
   endif ()


   if (ENABLE_CUDA)
      #
      # Retrieve compile flags for the current configuration
      # I haven't find a way to set host compiler flags for all the
      # possible configurations.
      #
      get_target_property( _amrex_flags_1 amrex COMPILE_OPTIONS)

      if (NOT CMAKE_CXX_FLAGS)
         get_target_property( _amrex_flags_2 Flags_CXX INTERFACE_COMPILE_OPTIONS)
      endif()

      set(_amrex_flags)
      if (_amrex_flags_1)
         list(APPEND _amrex_flags ${_amrex_flags_1})
      endif ()
      if (_amrex_flags_2)
         list(APPEND _amrex_flags ${_amrex_flags_2})
      endif ()

      eval_genex(_amrex_flags CXX  ${CMAKE_CXX_COMPILER_ID}
         COMP_VERSION ${CMAKE_CXX_COMPILER_VERSION}
         CONFIG       ${CMAKE_BUILD_TYPE}
         INTERFACE    BUILD
         STRING )

      if (_amrex_flags)
         target_compile_options(amrex PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${_amrex_flags}>)
      endif ()

      #
      # Add manually nvToolsExt if tiny profiler or base profiler are on.n
      # CMake >= 3.17 provides the module FindCUDAToolkit to do this natively.
      #
      if (ENABLE_TINY_PROFILE OR ENABLE_BASE_PROFILE )
          find_library(LIBNVTOOLSEXT nvToolsExt PATHS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
          target_link_libraries(amrex PUBLIC ${LIBNVTOOLSEXT})
      endif ()

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
        set(_${_ll}${_p} "${${_p}}")
        eval_genex( _${_ll}${_p} ${_l} ${CMAKE_${_l}_COMPILER_ID}
           COMP_VERSION ${CMAKE_${_l}_COMPILER_VERSION}
           CONFIG       ${CMAKE_BUILD_TYPE}
           INTERFACE    BUILD)

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
   if (CMAKE_Fortran_COMPILER_LOADED)
      message( STATUS "   Fortran compiler         = ${CMAKE_Fortran_COMPILER}")
   endif ()
   if (ENABLE_CUDA)
      message( STATUS "   CUDA compiler            = ${CMAKE_CUDA_COMPILER}")
   endif ()

   message( STATUS "   C++ defines              = ${_cxx_defines}")
   if (CMAKE_Fortran_COMPILER_LOADED)
      message( STATUS "   Fortran defines          = ${_fortran_defines}")
   endif ()

   message( STATUS "   C++ flags                = ${_cxx_flags}")
   if (CMAKE_Fortran_COMPILER_LOADED)
      message( STATUS "   Fortran flags            = ${_fortran_flags}")
   endif ()
   if (ENABLE_CUDA)
      message( STATUS "   CUDA flags               = ${CMAKE_CUDA_FLAGS_${AMREX_BUILD_TYPE}} ${CMAKE_CUDA_FLAGS}"
         "${AMREX_CUDA_FLAGS}")
   endif ()
   message( STATUS "   C++ include paths        = ${_cxx_includes}")
   if (CMAKE_Fortran_COMPILER_LOADED)
      message( STATUS "   Fortran include paths    = ${_fortran_includes}")
   endif ()
   message( STATUS "   Link line                = ${_cxx_link_line}")

endfunction ()
