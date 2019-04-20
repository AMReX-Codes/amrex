# 
#
# FUNCTION:  set_compiler_flags_preset
# 
# Set the compiler flags for target "_target".
# This function requires target "_target" to be already existent.
#
# Author: Michele Rosso
# Date  : Mar 6, 2019
#
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
   set(_cxx_gnu_dbg     "$<AND:${_cxx},${_debug},${_gnu}>")
   set(_cxx_gnu_rel     "$<AND:${_cxx},${_release},${_gnu}>")

   set(_intel           "$<CXX_COMPILER_ID:Intel>")
   set(_cxx_intel       "$<AND:${_cxx},${_intel}>")
   set(_cxx_intel_dbg   "$<AND:${_cxx},${_debug},${_intel}>")
   set(_cxx_intel_rel   "$<AND:${_cxx},${_release},${_intel}>")

   set(_pgi             "$<CXX_COMPILER_ID:PGI>")
   set(_cxx_pgi         "$<AND:${_cxx},${_pgi}>")
   set(_cxx_pgi_dbg     "$<AND:${_cxx},${_debug},${_pgi}>")
   set(_cxx_pgi_rel     "$<AND:${_cxx},${_release},${_pgi}>")

   set(_cray            "$<CXX_COMPILER_ID:Cray>")
   set(_cxx_cray        "$<AND:${_cxx},${_cray}>")
   set(_cxx_cray_dbg    "$<AND:${_cxx},${_debug},${_cray}>")
   set(_cxx_cray_rel    "$<AND:${_cxx},${_release},${_cray}>")

   set(_clang           "$<CXX_COMPILER_ID:Clang>")
   set(_cxx_clang       "$<AND:${_cxx},${_clang}>")
   set(_cxx_clang_dbg   "$<AND:${_cxx},${_debug},${_clang}>")
   set(_cxx_clang_rel   "$<AND:${_cxx},${_release},${_clang}>")
   
   set(_apple           "$<CXX_COMPILER_ID:AppleClang>")
   set(_cxx_apple       "$<AND:${_cxx},${_apple}>")
   set(_cxx_apple_dbg   "$<AND:${_cxx},${_debug},${_apple}>")
   set(_cxx_apple_rel   "$<AND:${_cxx},${_release},${_apple}>")

   
   # #
   # # Set REQUIRED fortran flags (required by at least amrdata)
   # # 
   # target_compile_options ( ${_target}
   #    PRIVATE
   #    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:$<${_fortran}:-ffixed-line-length-none -ffree-line-length-none>> 
   #    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:$<${_fortran}:-extend_source>>
   #    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Cray">:$<${_fortran}:-N 255 -h list=a>>
   #    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","PGI">:$<${_fortran}:-Mextend>> )


   # C++ REQUIRED flags
   # Until "cxx_std_11" and similar options are available (CMake >= 3.8 )
   # add c++11 support manually in order to have transitive property
   # Maybe here we should use the NOT genex 
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




