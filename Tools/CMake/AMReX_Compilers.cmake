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
   # Load defaults
   # 
   load_compiler_defaults ()

   # Overwrite defaults if user has set CMAKE_<lang>_FLAGS
   # Not the standard CMake-way of doing things
   if ( NOT CMAKE_Fortran_FLAGS )
      set ( CMAKE_Fortran_FLAGS ${AMREX_Fortran_FLAGS_${CMAKE_BUILD_TYPE}} )
   endif()

   if ( NOT CMAKE_CXX_FLAGS )
      set ( CMAKE_CXX_FLAGS ${AMREX_CXX_FLAGS_${CMAKE_BUILD_TYPE}} )
   endif() 

   target_compile_options ( amrex
      PUBLIC
      "$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_FLAGS}>" 
      "$<$<COMPILE_LANGUAGE:CXX>:${CMAKE_CXX_FLAGS}>" )
   
   if (DEFINED ENABLE_FPE)
      if (ENABLE_FPE)
	 target_compile_options ( amrex
	    PUBLIC
	    "$<$<COMPILE_LANGUAGE:Fortran>:${AMREX_Fortran_FLAGS_FPE}>" 
	    "$<$<COMPILE_LANGUAGE:CXX>:${AMREX_CXX_FLAGS_FPE}>" ) 
      endif ()
   else ()
      message (AUTHOR_WARNING "Variable ENABLE_FPE is not defined")
   endif ()
   
endfunction () 


# 
#
# MACRO:  load_compiler_defaults
# 
# Set the following variables to the default compiler flags:
#
# AMREX_<lang>_FLAGS_<CMAKE_BUILD_TYPE>
# AMREX_<lang>_FLAGS_FPE
# 
# where <lang> is either C or Fortran.
#
# The *_FPE flags enable floating point exceptions
#
# Author: Michele Rosso
# Date  : June 26, 2018
#
# 
macro ( load_compiler_defaults )

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
   # Load defaults
   # 
   if ( ${CMAKE_C_COMPILER_ID} STREQUAL "GNU")  ### GNU compiler ###
      
      set ( AMREX_Fortran_FLAGS_Release -O3)
      set ( AMREX_Fortran_FLAGS_Debug
	 -g -O0 -ggdb -fbounds-check -fbacktrace -Wuninitialized -Wunused -finit-real=snan -finit-integer=2147483647)
      set ( AMREX_Fortran_FLAGS_FPE -ffpe-trap=invalid,zero -ftrapv)
      

      set ( AMREX_CXX_FLAGS_Debug -g -O0 -fno-inline -ggdb -Wall -Wno-sign-compare)
      set ( AMREX_CXX_FLAGS_Release -O3 )
      set ( AMREX_CXX_FLAGS_FPE -ftrapv )
            
   elseif ( ${CMAKE_C_COMPILER_ID} STREQUAL "Intel")  ### Intel compiler ###

      set ( AMREX_Fortran_FLAGS_Debug -g -O0 -traceback -check bounds,uninit,pointers)
      set ( AMREX_Fortran_FLAGS_Release -O3 -ip -qopt-report=5 -qopt-report-phase=vec)
      set ( AMREX_Fortran_FLAGS_FPE )
      
      set ( AMREX_CXX_FLAGS_Debug     -g -O0 -traceback -Wcheck)
      set ( AMREX_CXX_FLAGS_Release   -O3 -ip -qopt-report=5 -qopt-report-phase=vec)
      set ( AMREX_CXX_FLAGS_FPE )
      
   elseif ( ${CMAKE_C_COMPILER_ID} STREQUAL "PGI")  ### PGI compiler ###

      set ( AMREX_Fortran_FLAGS_Debug -g -O0 -Mbounds -Ktrap=divz,inv -Mchkptr)
      set ( AMREX_Fortran_FLAGS_Release -gopt -fast)
      set ( AMREX_Fortran_FLAGS_FPE )
   
      set ( AMREX_CXX_FLAGS_Debug -O0 -Mbounds)
      set ( AMREX_CXX_FLAGS_Release -gopt -fast)
      set ( AMREX_CXX_FLAGS_FPE )

   elseif ( ${CMAKE_C_COMPILER_ID} STREQUAL "Cray")  ### Cray compiler ###

      set ( AMREX_Fortran_FLAGS_Debug -g -O0 -e i)
      set ( AMREX_Fortran_FLAGS_Release -O2)
      set ( AMREX_Fortran_FLAGS_FPE )
      
      set ( AMREX_CXX_FLAGS_Debug -g -O0)
      set ( AMREX_CXX_FLAGS_Release -O2)
      set ( AMREX_CXX_FLAGS_FPE )
      
   elseif ()

      set ( AMREX_Fortran_FLAGS_Debug )
      set ( AMREX_Fortran_FLAGS_Release )
      set ( AMREX_Fortran_FLAGS_FPE )
      
      set ( AMREX_CXX_FLAGS_Debug )
      set ( AMREX_CXX_FLAGS_Release )
      set ( AMREX_CXX_FLAGS_FPE )

      message ( WARNING "Compiler NOT recognized: ID is ${CMAKE_C_COMPILER_ID}. No default flags are loaded!" )
      
   endif ()
   
endmacro ()




