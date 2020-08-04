#
# This module defines INTERFACE targets defining presents of compiler flags:
# The INTERFACE targets defined here are:
#
#   Flags_CXX                 --> Optional flags for C++ code
#   Flags_Fortran             --> Optional flags for Fortran code
#   Flags_FPE                 --> Floating-Point Exception flags for both C++ and Fortran
#
# These INTERFACE targets can be added to the AMReX export set.
#
include_guard(GLOBAL)

#
# Create helper variables containing genex
# The variables created are:
#
#     _<lang>_<id>
#     _<lang>_<id>_dbg
#     _<lang>_<id>_rel
#
# for every combination of
#
#     <lang> = cxx,fortran
#     <id>   = gnu,intel,pgi,cray,clang,appleclang,msvc
#
foreach (_language CXX Fortran )
   set(_comp_lang   "$<COMPILE_LANGUAGE:${_language}>")
   string(TOLOWER "${_language}" _lang)

   foreach (_comp GNU Intel PGI Cray Clang AppleClang MSVC )
      string(TOLOWER "${_comp}" _id)
      # Define variables
      set(_comp_id              "$<${_language}_COMPILER_ID:${_comp}>")
      set(_${_lang}_${_id}      "$<AND:${_comp_lang},${_comp_id}>")
      set(_${_lang}_${_id}_dbg  "$<AND:${_comp_lang},${_comp_id},$<CONFIG:Debug>>")
      set(_${_lang}_${_id}_rel  "$<AND:${_comp_lang},${_comp_id},$<CONFIG:Release>>")
      unset(_comp_id)
   endforeach ()

   unset(_comp_lang)
   unset(_lang)
endforeach ()


#
# C++ flags
#
add_library(Flags_CXX INTERFACE)
add_library(AMReX::Flags_CXX ALIAS Flags_CXX)

target_compile_options( Flags_CXX
   INTERFACE
   $<${_cxx_gnu_dbg}:-O0 -ggdb -Wall -Wno-sign-compare -Wno-unused-but-set-variable -Werror=return-type>
#    $<$<VERSION_GREATER:$<CXX_COMPILER_VERSION>,5.0>:-Wnull-dereference>
   $<${_cxx_gnu_rel}:-Werror=return-type>
   $<${_cxx_intel_dbg}:-O0 -traceback -Wcheck>
   $<${_cxx_intel_rel}:-ip -qopt-report=5 -qopt-report-phase=vec>
   $<${_cxx_pgi_dbg}:-O0 -Mbounds>
   $<${_cxx_pgi_rel}:-gopt -fast>
   $<${_cxx_cray_dbg}:-O0>
   $<${_cxx_cray_rel}:>
   $<${_cxx_clang_dbg}:-O0 -Wall -Wextra -Wno-sign-compare -Wno-unused-parameter -Wno-unused-variable>
   $<${_cxx_clang_rel}:>
   $<${_cxx_appleclang_dbg}:-O0 -Wall -Wextra -Wno-sign-compare -Wno-unused-parameter -Wno-unused-variable>
   $<${_cxx_appleclang_rel}:>
   )

#
# Fortran flags
#
add_library(Flags_Fortran INTERFACE)
add_library(AMReX::Flags_Fortran ALIAS Flags_Fortran)

target_compile_options( Flags_Fortran
   INTERFACE
   $<${_fortran_gnu_dbg}:-O0 -ggdb -fcheck=bounds -fbacktrace -Wuninitialized -Wunused
   -finit-real=snan -finit-integer=2147483647 -fimplicit-none>
   $<${_fortran_gnu_rel}:-fimplicit-none>
   $<${_fortran_intel_dbg}:-O0 -traceback -check bounds,uninit,pointers -implicitnone>
   $<${_fortran_intel_rel}:-ip -qopt-report=5 -qopt-report-phase=vec -implicitnone>
   $<${_fortran_pgi_dbg}:-O0 -Mbounds>
   $<${_fortran_pgi_rel}:-gopt -fast>
   $<${_fortran_cray_dbg}:-O0 -e i>
   $<${_fortran_cray_rel}:>
   )


#
# Floating point exceptions
#
add_library(Flags_FPE INTERFACE)
add_library(AMReX::Flags_FPE ALIAS Flags_FPE)

target_compile_options ( Flags_FPE
   INTERFACE
   $<${_fortran_gnu}:-ffpe-trap=invalid,zero -ftrapv>
   $<${_cxx_gnu}:-ftrapv>
   $<${_fortran_intel}:-fpe3>
   $<${_cxx_intel}:-fpe3>
   $<${_fortran_pgi}:-Ktrap=divz,inv>
   $<${_cxx_pgi}:>
   $<${_fortran_cray}:-K trap=fp>
   $<${_cxx_cray}:-K trap=fp>
   $<${_fortran_clang}:>
   $<${_cxx_clang}:-ftrapv>
   )

#
# Unset all the variables defined in this module
#
foreach (_lang cxx fortran)
   foreach (_comp gnu intel pgi cray clang appleclang msvc)
      unset(_${_lang}_${_comp})
      unset(_${_lang}_${_comp}_dbg)
      unset(_${_lang}_${_comp}_rel)
   endforeach ()
endforeach ()
