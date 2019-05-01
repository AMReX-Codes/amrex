#
# This module defines INTERFACE targets defining presents of compiler flags:
# The INTERFACE targets defined here are:
#
#   Flags_CXX                 --> Optional flags for C++ code
#   Flags_Fortran             --> Optional flags for Fortran code
#   Flags_Fortran_REQUIRED    --> Required Fortran flags for some components of AMReX
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
#     <id>   = gnu,intel,pgi,cray,clang,appleclang 
#
foreach (_language CXX Fortran )
   set(_comp_lang   "$<COMPILE_LANGUAGE:${_language}>")
   string(TOLOWER "${_language}" _lang)   

   foreach (_comp GNU Intel PGI Cray Clang AppleClang )
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

target_compile_options( Flags_CXX
   INTERFACE
   $<${_cxx_gnu_dbg}:-O0 -fno-inline -ggdb -Wall -Wno-sign-compare>
   $<${_cxx_gnu_rel}:>
   $<${_cxx_intel_dbg}:-O0 -traceback -Wcheck>
   $<${_cxx_intel_rel}:-ip -qopt-report=5 -qopt-report-phase=vec>
   $<${_cxx_pgi_dbg}:-O0 -Mbounds>
   $<${_cxx_pgi_rel}:-gopt -fast>
   $<${_cxx_cray_dbg}:-O0>
   $<${_cxx_cray_rel}:>
   )

#
# Fortran flags
# 
add_library(Flags_Fortran INTERFACE)

target_compile_options( Flags_Fortran
   INTERFACE
   $<${_fortran_gnu_dbg}:-O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized -Wunused
   -finit-real=snan -finit-integer=2147483647>
   $<${_fortran_gnu_rel}:>
   $<${_fortran_intel_dbg}:-O0 -traceback -check bounds,uninit,pointers>
   $<${_fortran_intel_rel}:-ip -qopt-report=5 -qopt-report-phase=vec>
   $<${_fortran_pgi_dbg}:-O0 -Mbounds>
   $<${_fortran_pgi_rel}:-gopt -fast>
   $<${_fortran_cray_dbg}:-O0 -e i>
   $<${_fortran_cray_rel}:>
   )

#
# Fortran REQUIRED flags -- This is for internal use only: it useless to export it
#
add_library(Flags_Fortran_REQUIRED INTERFACE)

target_compile_options( Flags_Fortran_REQUIRED
   INTERFACE
   $<${_fortran_gnu}:-ffixed-line-length-none -ffree-line-length-none>
   $<${_fortran_intel}:-extend_source>
   $<${_fortran_pgi}:-N 255 -h list=a>
   $<${_fortran_cray}:-Mextend>
   )


#
# Floating point exceptions
# 
add_library(Flags_FPE INTERFACE)

target_compile_options ( Flags_FPE
   INTERFACE
   $<${_fortran_gnu}:-ffpe-trap=invalid,zero -ftrapv>
   $<${_cxx_gnu}:-ftrapv>
   $<${_fortran_intel}:-fpe3>
   $<${_cxx_intel}:-fpe3>
   $<${_fortran_cray}:-K trap=fp>
   $<${_cxx_cray}:-K trap=fp>
   $<${_fortran_pgi}:-Ktrap=divz,inv>
   $<${_cxx_pgi}:>
   )

#
# Unset all the variables defined in this module
#
foreach (_lang cxx fortran)
   foreach (_comp gnu intel pgi cray clang apple)
      unset(_${_lang}_${_comp})
      unset(_${_lang}_${_comp}_dbg)
      unset(_${_lang}_${_comp}_rel)     
   endforeach ()
endforeach ()

