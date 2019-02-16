# 
#
# FUNCTION: set_amrex_compiler_flags
# 
# Set the compiler flags for target "amrex".
# This function requires target "amrex" to be already existent.
#
# Author: Michele Rosso
# Date  : June 26, 2018
#
# 
function ( set_amrex_compiler_flags )

   # 
   # Check if target "amrex" has been defined before
   # calling this macro
   #
   if ( NOT TARGET amrex )
      message (FATAL_ERROR "Target 'amrex' must be defined before calling function 'set_amrex_compilers'" )
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
   # Set Fortran Flags only if not provided by user
   # Since genex "Fortran_COMPILER_ID" is not supported, we use genex
   # "STREQUAL" as workaround 
   #
   if ( NOT CMAKE_Fortran_FLAGS )     
       target_compile_options ( amrex
         PUBLIC
         # GNU Debug
         $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:$<$<COMPILE_LANGUAGE:Fortran>:
         -O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized -Wunused -finit-real=snan -finit-integer=2147483647>>>>
         # GNU Release
         $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:$<$<COMPILE_LANGUAGE:Fortran>:
         >>>>
         # Intel Debug
         $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:$<$<COMPILE_LANGUAGE:Fortran>:
         -O0 -traceback -check bounds,uninit,pointers>>>>
         # Intel Release
         $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:$<$<COMPILE_LANGUAGE:Fortran>:
         -ip -qopt-report=5 -qopt-report-phase=vec>>>>
         # Cray Debug
         $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Cray">:$<$<COMPILE_LANGUAGE:Fortran>:
         -O0 -e i>>>>
         # Cray Release 
         $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Cray">:$<$<COMPILE_LANGUAGE:Fortran>:
         >>>>
         # PGI Debug
         $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","PGI">:$<$<COMPILE_LANGUAGE:Fortran>:
         -O0 -Mbounds -Mchkptr>>>>
         # PGI Release
         $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","PGI">:$<$<COMPILE_LANGUAGE:Fortran>:
         -gopt -fast>>>>
         )            
   endif ()

   #
   # Set REQUIRED fortran flags 
   # 
   target_compile_options ( amrex
      PRIVATE
      # GNU 
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:$<$<COMPILE_LANGUAGE:Fortran>:
      -ffixed-line-length-none -ffree-line-length-none >> 
      # Intel 
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:$<$<COMPILE_LANGUAGE:Fortran>:
      -extend_source>>
      # Cray 
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Cray">:$<$<COMPILE_LANGUAGE:Fortran>:
      -N 255 -h list=a>>
      # PGI 
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","PGI">:$<$<COMPILE_LANGUAGE:Fortran>:
      -Mextend>> )

   #
   # Set C++ Flags only if not provided by user
   # 
   if ( NOT CMAKE_CXX_FLAGS )
      target_compile_options ( amrex
	 PUBLIC
	 # GNU Debug
	 $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<CXX_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	 -O0 -fno-inline -ggdb -Wall -Wno-sign-compare>>>>
	 # GNU Release
	 $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<CXX_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	 >>>>
	 # Intel Debug
	 $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<CXX_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:CXX>:
	 -O0 -traceback -Wcheck>>>>
	 # Intel Release
	 $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<CXX_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:CXX>:
	 -ip -qopt-report=5 -qopt-report-phase=vec>>>>
	 # Cray Debug
	 $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<CXX_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:CXX>:
	 -O0>>>>
	 # Cray Release 
	 $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<CXX_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:CXX>:
	 >>>>
	 # PGI Debug
	 $<BUILD_INTERFACE:$<$<CONFIG:Debug>:$<$<CXX_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:CXX>:
	 -O0 -Mbounds>>>>
	 # PGI Release
	 $<BUILD_INTERFACE:$<$<CONFIG:Release>:$<$<CXX_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:CXX>:
	 -gopt -fast>>>>
	 )	  
   endif ()


   # C++ REQUIRED flags
   # Until "cxx_std_11" and similar options are available (CMake >= 3.8 )
   # add c++11 support manually in order to have transitive property
   if (NOT ENABLE_3D_NODAL_MLMG)
      target_compile_options ( amrex
         PUBLIC
         $<$<CXX_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:CXX>:-h std=c++11 -h list=a>>
         $<$<CXX_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>>
         $<$<CXX_COMPILER_ID:Clang>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>>
         $<$<CXX_COMPILER_ID:AppleClang>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>>
         $<$<CXX_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>>
         $<$<CXX_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>> )
   else ()
      target_compile_options ( amrex
         PUBLIC
         $<$<CXX_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:CXX>:-h std=c++14 -h list=a>>
         $<$<CXX_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++14>>
         $<$<CXX_COMPILER_ID:Clang>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++14>>
         $<$<CXX_COMPILER_ID:AppleClang>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++14>>
         $<$<CXX_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++14>>
         $<$<CXX_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++14>> )
   endif ()

   
   #
   # Floating-point exceptions flags only if enabled
   # (I think these flags could be added tp both Fortran and
   # c++ without differentiating by language)
   # 
   if (DEFINED ENABLE_FPE)
      if (ENABLE_FPE)
   	 target_compile_options ( amrex
   	    PUBLIC
	    # GNU
	    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:$<$<COMPILE_LANGUAGE:Fortran>:
	    -ffpe-trap=invalid,zero -ftrapv>>
	    $<$<CXX_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	    -ftrapv>>
	    # Intel
	    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:$<$<COMPILE_LANGUAGE:Fortran>:
	    -fpe3>>
	    $<$<CXX_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:CXX>:
	    -fpe3>>
	    # Cray
	    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Cray">:$<$<COMPILE_LANGUAGE:Fortran>:
	    -K trap=fp>>
	    $<$<CXX_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:CXX>:
	    -K trap=fp>>
	    #  PGI
	    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","PGI">:$<$<COMPILE_LANGUAGE:Fortran>:
	    -Ktrap=divz,inv>>
	    $<$<CXX_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:CXX>:
	    >> )
      endif ()
   else ()
      message (AUTHOR_WARNING "Variable ENABLE_FPE is not defined")
   endif ()
   
endfunction () 




