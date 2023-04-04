#
#
# FUNCTION: get_amrex_version
#
# Retrieves AMReX version info and sets internal cache variables
# AMREX_GIT_VERSION, AMREX_RELEASE_NUMBER, and AMREX_PKG_VERSION
#
#
function (get_amrex_version)

   find_package(Git QUIET)

   set( _tmp "" )

   # Try to inquire software version from git
   if ( EXISTS ${CMAKE_CURRENT_LIST_DIR}/.git AND ${GIT_FOUND} )
      execute_process ( COMMAND git describe --abbrev=12 --dirty --always --tags
         WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
         OUTPUT_VARIABLE _tmp )
      string( STRIP "${_tmp}" _tmp )
      # filter invalid descriptions in shallow git clones
      if (NOT _tmp MATCHES "^([0-9]+)\\.([0-9]+)(\\.([0-9]+))*(-.*)*$")
         set( _tmp "")
      endif ()
   endif()

   # Grep first line from file CHANGES if cannot find version from Git
   if (NOT _tmp)
      file(STRINGS ${CMAKE_CURRENT_LIST_DIR}/CHANGES ALL_VERSIONS REGEX "#")
      list(GET ALL_VERSIONS 0 _tmp)
      string(REPLACE "#" "" _tmp "${_tmp}")
      string(STRIP "${_tmp}" _tmp )
      string(SUBSTRING "{_tmp}" 5 -1 _tmp2)
      if (NOT _tmp2) # don't want to add `.0` if it is already something like 21.10.2
         set(_tmp "${_tmp}.0")
      endif ()
      unset(_tmp2)
   endif ()

   set( AMREX_GIT_VERSION "${_tmp}" CACHE INTERNAL "" )
   unset(_tmp)

   # Package version is a modified form of AMREX_GIT_VERSION
   if (AMREX_GIT_VERSION)
      string(FIND "${AMREX_GIT_VERSION}" "-" _idx REVERSE)
      string(SUBSTRING "${AMREX_GIT_VERSION}" 0 "${_idx}" _pkg_version )
      string(FIND "${_pkg_version}" "-" _idx REVERSE)
      string(SUBSTRING "${_pkg_version}" 0 "${_idx}" _pkg_version )
      string(REPLACE "-" "." _pkg_version "${_pkg_version}")
   endif ()

   set( AMREX_PKG_VERSION "${_pkg_version}" CACHE INTERNAL "" )
   unset(_pkg_version)

   if (AMREX_GIT_VERSION)
      string(FIND "${AMREX_GIT_VERSION}" "-" _idx)
      string(SUBSTRING "${AMREX_GIT_VERSION}" 0 "${_idx}" _rel_number )
      string(REPLACE "." "" _rel_number "${_rel_number}")
      string(SUBSTRING "${_rel_number}" 0 4 _rel_yymm)
      string(SUBSTRING "${_rel_number}" 4 2 _rel_patch)
      string(LENGTH "${_rel_patch}" _rel_patch_len)
      if (_rel_patch_len EQUAL 0)
         string(PREPEND _rel_patch "00")
      elseif (_rel_patch_len EQUAL 1)
         string(PREPEND _rel_patch "00")
      endif ()
      string(CONCAT _rel_number "${_rel_yymm}" "${_rel_patch}")
      unset(_rel_yymm)
      unset(_rel_patch)
      unset(_rel_patch_len)
   endif ()

   set( AMREX_RELEASE_NUMBER "${_rel_number}" CACHE INTERNAL "" )
   unset(_rel_number)

endfunction ()


#
#
# FUNCTION: print
#
# Debug function to print a variable
#
# Arguments:
#
#    _var = the variable to print
#
#
function ( print _var )
   message(" ${_var} = ${${_var}}" )
endfunction ()

#
#
# FUNCTION: print_list
#
# Debug function to print a list as an columns of entries
#
# Arguments:
#
#    _list = the list to print
#
#
function ( print_list _list )

   list( LENGTH ${_list} _len )

   if ( ${_len} GREATER 0 )
      message("")
      message( STATUS "LIST NAME:  ${_list}")
      foreach ( _item ${${_list}})
         message ( STATUS "  ${_item}")
      endforeach ()
      message("")
   endif ()

endfunction ()


#
#
# MACRO: set_default_config_flags
#
#
# Set CMake_<LANG>_FLAGS_<CONFIG> to default values
# if the variable is empty
#
macro ( set_default_config_flags )

  if ( NOT CMAKE_Fortran_FLAGS_DEBUG )
    set (CMAKE_Fortran_FLAGS_DEBUG "-g")
  endif ()

  if ( NOT CMAKE_Fortran_FLAGS_RELEASE )
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  endif ()

  if ( NOT CMAKE_CXX_FLAGS_DEBUG )
    set (CMAKE_CXX_FLAGS_DEBUG "-g")
  endif ()

  if ( NOT CMAKE_CXX_FLAGS_RELEASE )
    set (CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
  endif ()

endmacro ()


#
#
# FUNCTION: set_mininum_compiler_version
#
# Check whether the compiler version for a certain lnaguage is
# >= a required minimum. If not, stop with a fatal error
#
# Arguments:
#
#    _lang            = the language
#    _comp_id         = the compiler ID
#    _minimum_version = the minimum version required for _comp_id
#
#
function (set_mininum_compiler_version _lang _comp_id  _minimum_version)
   if (  (CMAKE_${_lang}_COMPILER_ID STREQUAL _comp_id ) AND
         (CMAKE_${_lang}_COMPILER_VERSION VERSION_LESS _minimum_version ) )
      message( FATAL_ERROR
         "\n${_comp_id} compiler version is ${CMAKE_${_lang}_COMPILER_VERSION}."
         " Minimum required is ${_minimum_version}.\n")
   endif ()
endfunction ()



#
#
# FUNCTION: check_cuda_host_compiler
#
#
# Makes sure the CUDA host compiler and CXX compiler are the same.
# CMake lets you decide which host compiler to use via the env variable
# CUDAHOSTCXX and the CMake variable CMAKE_CUDA_HOST_COMPILER.
# For the time being we force the CUDA host compiler to be the C++ compiler.
#
# Note: just comparing the CMAKE_..._COMPILER vars is not sufficient and raises
#       false negatives on e.g. /usr/bin/g++-8 and /usr/bin/c++
# Note: blocked by https://gitlab.kitware.com/cmake/cmake/-/issues/20901
#
#
function (check_cuda_host_compiler)
   if (CMAKE_CUDA_HOST_COMPILER)
      if (NOT "${CMAKE_CXX_COMPILER}" STREQUAL "${CMAKE_CUDA_HOST_COMPILER}")
         if (NOT "$ENV{CUDAHOSTCXX}" STREQUAL "" OR NOT "$ENV{CXX}" STREQUAL "")
            message(WARNING "CUDA host compiler "
                            "(${CMAKE_CUDA_HOST_COMPILER}) "
                            "does not match the C++ compiler "
                            "(${CMAKE_CXX_COMPILER})! "
                            "Consider setting the CXX and CUDAHOSTCXX environment "
                            "variables.")
         endif ()
      endif ()
   endif ()
endfunction ()


#
#
# FUNCTION: convert_cuda_archs
#
#
# cuda_select_nvcc_arch_flags accepts CUDA architecture in the form of
# names (Turing, Volta, ...) or decimal numbers (10.0, 9.0, ...).
# However, CMAKE_CUDA_ARCHITECTURES only accepts integer numbers.
# We need to convert the latter to decimal format, else cuda_select_nvcc_arch_flags
# will complain
#
# Arguments:
#
#    _cuda_archs   = the target architecture(s)
#
#
function (convert_cuda_archs  _cuda_archs)

   foreach (_item IN LISTS ${_cuda_archs})
      string(REGEX MATCH "\\." _has_decimal "${_item}")
      string(REGEX MATCH "[0-9]+" _is_number "${_item}")

      if (NOT _has_decimal AND _is_number)
         math(EXPR _int "${_item}/10" OUTPUT_FORMAT DECIMAL)
         math(EXPR _mod "${_item}%10" OUTPUT_FORMAT DECIMAL)
         list(APPEND _tmp "${_int}.${_mod}")
      else ()
         list(APPEND _tmp ${_item})
      endif()
   endforeach ()

   set(${_cuda_archs} ${_tmp} PARENT_SCOPE)

endfunction ()


#
#
# FUNCTION: set_cuda_architectures
#
#
# Detects the cuda capabilities of the GPU and set the internal
# variable AMREX_CUDA_ARCHS.
#
# Arguments:
#
#    _cuda_archs   = the target architecture(s) (select "Auto" for autodetection)
#
# Note: if no target arch is specified, it will try to determine
# which GPU architecture is supported by the system. If more than one is found,
# it will build for all of them.
# If autodetection fails, a list of "common" architectures is assumed.
#
function (set_cuda_architectures _cuda_archs)

   set(_archs ${${_cuda_archs}})
   convert_cuda_archs(_archs)

   include(FindCUDA/select_compute_arch)
   cuda_select_nvcc_arch_flags(_nvcc_arch_flags ${_archs})

   # Extract architecture number: anything less than 6.0 must go
   string(REPLACE "-gencode;" "-gencode=" _nvcc_arch_flags "${_nvcc_arch_flags}")

   foreach (_item IN LISTS _nvcc_arch_flags)
      # Match one time the regex [0-9]+.
      # [0-9]+ means any number between 0 and 9 will be matched one or more times (option +)
      string(REGEX MATCH "[0-9]+" _cuda_compute_capability "${_item}")

      if (_cuda_compute_capability LESS 60)
         message(STATUS "Ignoring unsupported CUDA architecture ${_cuda_compute_capability}")
      else ()
         list(APPEND _tmp ${_cuda_compute_capability})
      endif ()
   endforeach ()

   set(AMREX_CUDA_ARCHS ${_tmp} CACHE INTERNAL "CUDA archs AMReX is built for")

endfunction()


#
#
# FUNCTION: set_nvcc_arch_flags
#
#
# Detects the cuda capabilities of the GPU and set the internal
# variable NVCC_ARCH_FLAGS.
#
# Arguments:
#
#    _cuda_archs   = the target architecture(s) (select "Auto" for autodetection)
#    _lto          = true if LTO flags are required
#
# Note: if no target arch is specified, it will try to determine
# which GPU architecture is supported by the system. If more than one is found,
# it will build for all of them.
# If autodetection fails, a list of “common” architectures is assumed.
#
function (set_nvcc_arch_flags _cuda_archs _lto)

   set(_archs ${${_cuda_archs}})
   convert_cuda_archs(_archs)

   include(FindCUDA/select_compute_arch)
   cuda_select_nvcc_arch_flags(_nvcc_arch_flags ${_archs})

   #
   # Remove unsupported architecture: anything less than 3.5 must go
   #
   string(REPLACE "-gencode;" "-gencode=" _nvcc_arch_flags "${_nvcc_arch_flags}")

   foreach (_item IN LISTS _nvcc_arch_flags)
      # Match one time the regex [0-9]+.
      # [0-9]+ means any number between 0 and 9 will be matched one or more times (option +)
      string(REGEX MATCH "[0-9]+" _cuda_compute_capability "${_item}")

      if (_cuda_compute_capability LESS 35)
         message(STATUS "Ignoring unsupported CUDA architecture ${_cuda_compute_capability}")
         list(REMOVE_ITEM _nvcc_arch_flags ${_item})
      endif ()

   endforeach ()

   if (${_lto})
      # we replace
      #   -gencode=arch=compute_NN,code=sm_NN
      # with
      #   -gencode=arch=compute_NN,code=lto_NN
      set(_nvcc_arch_flags_org ${_nvcc_arch_flags})
      foreach (_item IN LISTS _nvcc_arch_flags_org)
         string(REGEX MATCH "[0-9]+" _cuda_compute_capability "${_item}")
         string(REPLACE "code=sm_${_cuda_compute_capability}"
            "code=lto_${_cuda_compute_capability}"
            _nvcc_arch_flags "${_nvcc_arch_flags}")
      endforeach ()
   endif ()

   if (NOT _nvcc_arch_flags)
      message(FATAL_ERROR "the given target CUDA architectures are not supported by AMReX")
   endif ()

   #
   # Set architecture-dependent flags
   #
   string(REPLACE ";" " " _nvcc_arch_flags "${_nvcc_arch_flags}")
   set(NVCC_ARCH_FLAGS ${_nvcc_arch_flags} CACHE INTERNAL "CUDA architecture-dependent flags")

endfunction()
