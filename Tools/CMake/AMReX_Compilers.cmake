# 
#
# FUNCTION: set_amrex_compilers
# 
# Set the compiler flags for target "amrex".
# This function requires target "amrex" to be already existent.
#
# Author: Michele Rosso
# Date  : June 26, 2018
#
# 
function ( set_amrex_compilers )

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
   # Check the same compiler suite is used for all languages
   # 
   if (  NOT (${CMAKE_C_COMPILER_ID} STREQUAL ${CMAKE_Fortran_COMPILER_ID}) OR
	 NOT (${CMAKE_C_COMPILER_ID} STREQUAL ${CMAKE_CXX_COMPILER_ID}) )
      message ( FATAL_ERROR "C compiler ID does not match Fortran/C++ compiler ID" )
   endif ()

   
   #
   # Set Fortran Flags
   # 
   if ( CMAKE_Fortran_FLAGS )
      target_compile_options ( amrex
	 PUBLIC
	 $<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_FLAGS}> )
   else () # Set defaults
      target_compile_options ( amrex
	 PUBLIC
	 # GNU Debug
	 $<$<CONFIG:Debug>:$<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -g -O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized -Wunused -finit-real=snan -finit-integer=2147483647>>>
	 # GNU Release
	 $<$<CONFIG:Release>:$<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -O3>>>
	 # Intel Debug
	 $<$<CONFIG:Debug>:$<$<C_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -g -O0 -traceback -check bounds,uninit,pointers>>>
	 # Intel Release
	 $<$<CONFIG:Release>:$<$<C_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -O3 -ip -qopt-report=5 -qopt-report-phase=vec>>>
	 # Cray Debug
	 $<$<CONFIG:Debug>:$<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -g -O0 -e i >>>
	 # Cray Release 
	 $<$<CONFIG:Release>:$<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -O2>>>
	 # PGI Debug
	 $<$<CONFIG:Debug>:$<$<C_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -g -O0 -Mbounds -Mchkptr>>>
	 # PGI Release
	 $<$<CONFIG:Release>:$<$<C_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:Fortran>:
	 -gopt -fast>>>
	 )	  
   endif ()

   #
   # Set REQUIRED fortran flags
   # 
   target_compile_options ( amrex
      PRIVATE
      # GNU 
      $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
      -ffixed-line-length-none -ffree-line-length-none -fno-range-check -fno-second-underscore>>
      # Intel 
      $<$<C_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:Fortran>:
      -extend_source>>
      # Cray 
      $<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:Fortran>:
      -N 255 -h list=a>>
      # PGI 
      $<$<C_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:Fortran>:
      -extend>> )

   #
   # Set C++ Flags
   # 
   if ( CMAKE_CXX_FLAGS )
      target_compile_options ( amrex
	 PUBLIC
	 $<$<COMPILE_LANGUAGE:CXX>:${CMAKE_CXX_FLAGS}> )
   else () # Set defaults
      target_compile_options ( amrex
	 PUBLIC
	 # GNU Debug
	 $<$<CONFIG:Debug>:$<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	 -g -O0 -fno-inline -ggdb -Wall -Wno-sign-compare>>>
	 # GNU Release
	 $<$<CONFIG:Release>:$<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	 -O3>>>
	 # Intel Debug
	 $<$<CONFIG:Debug>:$<$<C_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:CXX>:
	 -g -O0 -traceback -Wcheck>>>
	 # Intel Release
	 $<$<CONFIG:Release>:$<$<C_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:CXX>:
	 -O3 -ip -qopt-report=5 -qopt-report-phase=vec>>>
	 # Cray Debug
	 $<$<CONFIG:Debug>:$<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:CXX>:
	 -g -O0>>>
	 # Cray Release 
	 $<$<CONFIG:Release>:$<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:CXX>:
	 -O2>>>
	 # PGI Debug
	 $<$<CONFIG:Debug>:$<$<C_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:CXX>:
	 -O0 -Mbounds>>>
	 # PGI Release
	 $<$<CONFIG:Release>:$<$<C_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:CXX>:
	 -gopt -fast>>>
	 )	  
   endif ()


   # C++ REQUIRED flags
   # Until "cxx_std_11" and similar options are available (CMake >= 3.8 )
   # add c++11 support manually in order to have transitive property
   target_compile_options ( amrex PUBLIC
      $<$<C_COMPILER_ID:Cray>:$<$<COMPILE_LANGUAGE:CXX>:-h std=c++11 -h list=a>>
      $<$<C_COMPILER_ID:PGI>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>>
      $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>>
      $<$<C_COMPILER_ID:Intel>:$<$<COMPILE_LANGUAGE:CXX>:-std=c++11>> )

   
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
	    $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
	    -ffpe-trap=invalid,zero -ftrapv>>
	    $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	    -ftrapv>>
	    # Intel
	    $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
	    -fpe3>>
	    $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	    -fpe3>>
	    # Cray
	    $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
	    -K trap=fp>>
	    $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	    -K trap=fp>>
	    #  PGI
	    $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:Fortran>:
	    -Ktrap=divz,inv>>
	    $<$<C_COMPILER_ID:GNU>:$<$<COMPILE_LANGUAGE:CXX>:
	    >> )
      endif ()
   else ()
      message (AUTHOR_WARNING "Variable ENABLE_FPE is not defined")
   endif ()
   
endfunction () 




