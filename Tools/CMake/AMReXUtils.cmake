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

   # Grep first line from file CHANGES.md if cannot find version from Git
   if (NOT _tmp)
      file(STRINGS ${CMAKE_CURRENT_LIST_DIR}/CHANGES.md ALL_VERSIONS REGEX "#")
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
# Expand the CUDA architecture aliases "all" and "all-major" into the concrete
# architectures they stand for, spelled the way the CUDA_ARCHITECTURES property does:
# SASS code for every architecture and PTX for the newest one only.
#
# The architectures a toolkit supports are queried from the CUDA compiler itself, so that
# the answer stays correct for a toolkit newer than the CMake release in use, whose
# CMAKE_CUDA_ARCHITECTURES_ALL[_MAJOR] would list architectures the compiler has dropped
# and miss the ones it added. Those variables are the fallback for a compiler that does not
# answer the query (e.g. clang -x cuda), which is also the list CMake itself expands the
# alias with in that case.
#
#   _alias  "all" or "all-major"
#   _var    name of a variable that receives the architecture list in the caller; empty if
#           neither the compiler nor CMake knows the architectures
#
function (amrex_expand_cuda_archs_alias _alias _var)
   set(_archs)
   execute_process(COMMAND ${CMAKE_CUDA_COMPILER} --list-gpu-arch
      OUTPUT_VARIABLE _out RESULT_VARIABLE _rv
      ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
   if (_rv EQUAL 0)
      string(REGEX MATCHALL "compute_([0-9]+)" _archs "${_out}")
      list(TRANSFORM _archs REPLACE "compute_" "")
      # they are not printed in ascending order (e.g. 110 before 103 with CUDA 13.3), and
      # both the earliest and the newest one are needed below
      list(SORT _archs COMPARE NATURAL)
   endif ()

   if (_archs AND _alias STREQUAL "all-major")
      # "all-major" covers every major architecture (minor revision 0) plus the earliest
      # architecture the toolkit supports, which need not be a major one (sm_75 with CUDA 13)
      list(GET _archs 0 _earliest)
      set(_major "${_earliest}")
      foreach (_arch IN LISTS _archs)
         math(EXPR _minor "${_arch} % 10")
         if (_minor EQUAL 0 AND NOT _arch STREQUAL _earliest)
            list(APPEND _major "${_arch}")
         endif ()
      endforeach ()
      set(_archs "${_major}")
   endif ()

   if (_archs)
      list(POP_BACK _archs _newest)
      list(TRANSFORM _archs APPEND "-real")
      list(APPEND _archs "${_newest}")
   elseif (_alias STREQUAL "all")
      set(_archs "${CMAKE_CUDA_ARCHITECTURES_ALL}")
   else ()
      set(_archs "${CMAKE_CUDA_ARCHITECTURES_ALL_MAJOR}")
   endif ()

   set(${_var} "${_archs}" PARENT_SCOPE)
endfunction ()


#
# Drop architectures below the minimum compute capability AMReX supports, as the deprecated
# FindCUDA-based selection used to do. Without this an unsupported architecture fails much
# later: in CMake's CUDA compiler detection, or deep inside the device code compilation of
# the first source that uses e.g. atomicAdd(double*).
#
#   _var         name of a variable holding the architecture list; updated in the caller
#   _from_alias  TRUE if the list came from expanding native/all/all-major, where dropping
#                unsupported entries is expected and only worth a status message
#
function (amrex_filter_cuda_archs _var _from_alias)
   set(_ok)
   set(_low)
   foreach (_arch IN LISTS ${_var})
      set(_supported TRUE)
      if (_arch MATCHES "^([0-9]+)")
         if (CMAKE_MATCH_1 LESS 60)
            set(_supported FALSE)
         endif ()
      endif ()
      if (_supported)
         list(APPEND _ok "${_arch}")
      else ()
         list(APPEND _low "${_arch}")
      endif ()
   endforeach ()

   if (_low)
      if (NOT _ok)
         message(FATAL_ERROR
            "The requested CUDA architecture(s) ${_low} are not supported by AMReX, which "
            "requires compute capability 6.0 or higher. Set CMAKE_CUDA_ARCHITECTURES to a "
            "supported architecture, e.g. -DCMAKE_CUDA_ARCHITECTURES=80.")
      elseif (_from_alias)
         message(STATUS "   ignoring CUDA architectures below 6.0: ${_low}")
      else ()
         message(WARNING
            "Ignoring the requested CUDA architecture(s) ${_low}: AMReX requires compute "
            "capability 6.0 or higher. Building for ${_ok}.")
      endif ()
      set(${_var} "${_ok}" PARENT_SCOPE)
   endif ()
endfunction ()
