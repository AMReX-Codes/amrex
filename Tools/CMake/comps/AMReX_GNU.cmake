#
# Setup for GNU compilers suite.
# This file should be included by 'setup_amrex_compilers()' only
# It assume that target 'amrex' and option 'ENABLE_FPE' have been
# already defined as well as 
#

# C++
set(cxx_gnu      "$<AND:$<CXX_COMPILER_ID:GNU>,$<COMPILE_LANGUAGE:CXX>>")
set(cxx_gnu_dbg  "$<AND:$<CONFIG:Debug>,${cxx_gnu}>")
set(cxx_gnu_rel  "$<AND:$<CONFIG:Release>,${cxx_gnu}>")

# Fortran
set(fortran_gnu      "$<AND:$<STREQUAL:\"${CMAKE_Fortran_COMPILER_ID}\",\"GNU\">,$<COMPILE_LANGUAGE:Fortran>>" )
set(fortran_gnu_dbg  "$<AND:$<CONFIG:Debug>,${fortran_gnu}>")
set(fortran_gnu_rel  "$<AND:$<CONFIG:Release>,${fortran_gnu}>")

#
# Set Fortran Flags only if not provided by user
#
if ( NOT CMAKE_Fortran_FLAGS )     
   target_compile_options ( amrex
      PUBLIC
      $<${fortran_gnu_dbg}:$<BUILD_INTERFACE:-O0 -ggdb -fcheck=bounds -fbacktrace -Wuninitialized
      -Wunused -finit-real=snan -finit-integer=2147483647>>
      $<${fortran_gnu_rel}:$<BUILD_INTERFACE:>>
      )            
endif ()

#
# Set REQUIRED fortran flags 
# 
target_compile_options( amrex
   PRIVATE
   $<${fortran_gnu}:-ffixed-line-length-none -ffree-line-length-none> )

#
# Set C++ Flags only if not provided by user
# 
if ( NOT CMAKE_CXX_FLAGS )
   target_compile_options ( amrex
      PUBLIC
      $<${cxx_gnu_dbg}:$<BUILD_INTERFACE:-O0 -ggdb -Wshadow -Wall -Wno-sign-compare
      -Wno-unused-but-set-variable -Werror=return-type>>    
      $<${cxx_gnu_rel}:$<BUILD_INTERFACE:>>
      )	  
endif ()

#
# Floating-point exceptions flags only if enabled
# (I think these flags could be added to both Fortran and
# c++ without differentiating by language)
# 
if (ENABLE_FPE)
   target_compile_options( amrex
      PUBLIC
      $<${fortran_gnu}:-ffpe-trap=invalid,zero -finit-real=snan -finit-integer=2147483647 -ftrapv>
      $<${cxx_gnu}:-ftrapv>
      )
endif ()

#
# Set compile features -- We just need to set C++ standard
# This could be done by using target_compile_features because it propagates
# the standard c++11 flag to the nvcc compiler and we do not want that since
# nvcc does not support c++14 
#
if (NOT ENABLE_3D_NODAL_MLMG)
   target_compile_options( amrex PUBLIC  $<${cxx_gnu}:-std=c++11> )
else ()
   target_compile_options( amrex PUBLIC  $<${cxx_gnu}:-std=c++14> $<${cuda_gnu}:-Xcompiler=-std=c++14>)
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

   target_compile_definitions( amrex PUBLIC
      BL_GCC_VERSION=${CMAKE_CXX_COMPILER_VERSION}
      BL_GCC_MAJOR_VERSION=${GCC_VERSION_MAJOR}
      BL_GCC_MINOR_VERSION=${GCC_VERSION_MINOR}
      )
endif ()
