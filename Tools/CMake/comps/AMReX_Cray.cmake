#
# Setup for CRAY compilers suite.
# This file should be included by 'setup_amrex_compilers()' only
# It assume that target 'amrex' and option 'ENABLE_FPE' have been
# already defined as well as 
#

# C++
set(cxx_cray      "$<AND:$<CXX_COMPILER_ID:Cray>,$<COMPILE_LANGUAGE:CXX>>")
set(cxx_cray_dbg  "$<AND:$<CONFIG:Debug>,${cxx_cray}>")
set(cxx_cray_rel  "$<AND:$<CONFIG:Release>,${cxx_cray}>")

# Fortran
set(fortran_cray      "$<AND:$<STREQUAL:\"${CMAKE_Fortran_COMPILER_ID}\",\"Cray\">,$<COMPILE_LANGUAGE:Fortran>>" )
set(fortran_cray_dbg  "$<AND:$<CONFIG:Debug>,${fortran_cray}>")
set(fortran_cray_rel  "$<AND:$<CONFIG:Release>,${fortran_cray}>")

#
# Set Fortran Flags only if not provided by user
#
if ( NOT CMAKE_Fortran_FLAGS )     
   target_compile_options ( amrex
      PUBLIC
      $<${fortran_cray_dbg}:$<BUILD_INTERFACE:-O0 -e i>>
      $<${fortran_cray_rel}:$<BUILD_INTERFACE:>>
      )            
endif ()

#
# Set REQUIRED fortran flags 
# 
target_compile_options( amrex PRIVATE  $<${fortran_cray}:-N 255 -h list=a> )

#
# Set C++ Flags only if not provided by user
# 
if ( NOT CMAKE_CXX_FLAGS )
   target_compile_options ( amrex
      PUBLIC
      $<${cxx_cray_dbg}:$<BUILD_INTERFACE:-O0>>
      $<${cxx_cray_rel}:$<BUILD_INTERFACE:>>
      )	  
endif ()


if (ENABLE_FPE)
   target_compile_options ( amrex
      PUBLIC
      $<${fortran_cray}:-K trap=fp>
      $<${cxx_cray}:-K trap=fp>
      )
endif ()


#
# Set compile features -- We just need to set C++ standard
# This could be done by using target_compile_features because it propagates
# the standard c++11 flag to the nvcc compiler and we do not want that since
# nvcc does not support c++14 
#
if (NOT ENABLE_3D_NODAL_MLMG)
   target_compile_options( amrex PUBLIC $<${cxx_cray}:-h std=c++11 -h list=a> )
else ()
   target_compile_options( amrex PUBLIC $<${cxx_cray}:-h std=c++14 -h list=a> )
endif ()
