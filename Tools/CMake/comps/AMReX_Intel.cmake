#
# Setup for Intel compilers suite.
# This file should be included by 'setup_amrex_compilers()' only
# It assume that target 'amrex' and option 'ENABLE_FPE' have been
# already defined as well as 
#

# C++
set(cxx_intel      "$<AND:$<CXX_COMPILER_ID:Intel>,$<COMPILE_LANGUAGE:CXX>>")
set(cxx_intel_dbg  "$<AND:$<CONFIG:Debug>,${cxx_intel}>")
set(cxx_intel_rel  "$<AND:$<CONFIG:Release>,${cxx_intel}>")

# Fortran
set(fortran_intel      "$<AND:$<STREQUAL:\"${CMAKE_Fortran_COMPILER_ID}\",\"Intel\">,$<COMPILE_LANGUAGE:Fortran>>" )
set(fortran_intel_dbg  "$<AND:$<CONFIG:Debug>,${fortran_intel}>")
set(fortran_intel_rel  "$<AND:$<CONFIG:Release>,${fortran_intel}>")

#
# Set Fortran Flags only if not provided by user
#
if ( NOT CMAKE_Fortran_FLAGS )     
   target_compile_options ( amrex
      PUBLIC
      $<${fortran_intel_dbg}:$<BUILD_INTERFACE:-O0 -traceback -check bounds,uninit,pointers>>
      $<${fortran_intel_rel}:$<BUILD_INTERFACE:-ip -qopt-report=5 -qopt-report-phase=vec>>
      )            
endif ()

#
# Set REQUIRED fortran flags 
# 
target_compile_options( amrex
   PRIVATE
   $<${fortran_intel}:-extend_source>
   )

#
# Set C++ Flags only if not provided by user
# 
if ( NOT CMAKE_CXX_FLAGS )
   target_compile_options ( amrex
      PUBLIC
      $<${cxx_intel_dbg}:$<BUILD_INTERFACE:-O0 -traceback -Wcheck>>
      $<${cxx_intel_rel}:$<BUILD_INTERFACE:-ip -qopt-report=5 -qopt-report-phase=vec>>
      )	  
endif ()

#
# Floating-point exceptions flags only if enabled
# (I think these flags could be added tp both Fortran and
# c++ without differentiating by language)
# 
if (ENABLE_FPE)
   target_compile_options( amrex PUBLIC
      $<${fortran_intel}:-fpe3>
      $<${cxx_intel}:-fpe3>
      )
endif ()

