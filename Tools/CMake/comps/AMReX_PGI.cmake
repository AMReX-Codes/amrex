#
# Setup for PGI compilers suite.
# This file should be included by 'setup_amrex_compilers()' only
# It assume that target 'amrex' and option 'ENABLE_FPE' have been
# already defined as well as 
#

# C++
set(cxx_pgi      "$<AND:$<CXX_COMPILER_ID:PGI>,$<COMPILE_LANGUAGE:CXX>>")
set(cxx_pgi_dbg  "$<AND:$<CONFIG:Debug>,${cxx_pgi}>")
set(cxx_pgi_rel  "$<AND:$<CONFIG:Release>,${cxx_pgi}>")

# Fortran
set(fortran_pgi      "$<AND:$<STREQUAL:\"${CMAKE_Fortran_COMPILER_ID}\",\"PGI\">,$<COMPILE_LANGUAGE:Fortran>>" )
set(fortran_pgi_dbg  "$<AND:$<CONFIG:Debug>,${fortran_pgi}>")
set(fortran_pgi_rel  "$<AND:$<CONFIG:Release>,${fortran_pgi}>")


#
# Set Fortran Flags only if not provided by user
# (CMake adds -Mbounds to fortran dor debug config)
#
if ( NOT CMAKE_Fortran_FLAGS )     
   target_compile_options ( amrex
      PUBLIC
      $<${fortran_pgi_dbg}:$<BUILD_INTERFACE:-Mchkptr>>
      $<${fortran_pgi_rel}:$<BUILD_INTERFACE:>>
      )            
endif ()

#
# Set REQUIRED fortran flags 
# 
target_compile_options( amrex PRIVATE $<${fortran_pgi}:-Mextend -Mdclchk>  )

#
# Set C++ Flags only if not provided by user
# 
if ( NOT CMAKE_CXX_FLAGS )
   target_compile_options ( amrex
      PUBLIC
      $<${cxx_pgi_dbg}:$<BUILD_INTERFACE:-Mbounds>>
      $<${cxx_pgi_rel}:$<BUILD_INTERFACE:>>
      )	  
endif ()


#
# Floating-point exceptions flags only if enabled
# (I think these flags could be added tp both Fortran and
# c++ without differentiating by language)
# 
if (ENABLE_FPE)
   target_compile_options ( amrex
      PUBLIC
      $<${fortran_pgi}:-Ktrap=divz,inv>
      $<${cxx_pgi}:>> )
endif ()

#
# Set compile features -- We just need to set C++ standard
# This could be done by using target_compile_features because it propagates
# the standard c++11 flag to the nvcc compiler and we do not want that since
# nvcc does not support c++14 
#
if (NOT ENABLE_3D_NODAL_MLMG)
   target_compile_options( amrex PUBLIC  $<${cxx_pgi}:-std=c++11> )
else ()
   target_compile_options( amrex PUBLIC  $<${cxx_pgi}:-std=c++14> )
endif ()

#
# CUDA
# 
if (ENABLE_CUDA)
   target_compile_options( amrex
      PUBLIC
      $<${fortran_pgi}:-Mcuda=cc${CUDA_ARCH},ptxinfo,fastmath,charstring>
      $<${fortran_pgi_dbg}:-Mcuda=debug>
      $<${fortran_pgi_rel}:-Mcuda=lineinfo>
      $<${fortran_pgi}:CUDA_HOME=${CUDA_HOME}>
      )

   if (CUDA_MAXREGCOUNT)
      target_compile_options( amrex
         PUBLIC
         $<${fortran_pgi}:-Mcuda=maxregcount:${CUDA_MAXREGCOUNT}>
         )
   endif ()

   if (ENABLE_CUDA_FORTRAN)
      target_compile_definitions( amrex PUBLIC
         $<${fortran_pgi}:AMREX_USE_CUDA_FORTRAN> )
   endif ()   
endif ()


#
# OpenAcc
#
if (ENABLE_ACC)
  # target_compile_options( amrex PUBLIC  $<${cxx_pgi}:-noacc>  $<${fortran_pgi}:-noacc>)
else ()
   target_compile_options( amrex PUBLIC  $<${cxx_pgi}:-noacc> $<${fortran_pgi}:-noacc>)   
endif ()

#
# Other flags
#
if (NOT ENABLE_FORTRAN_INTERFACES)
   target_compile_options( amrex PUBLIC $<${fortran_pgi}:-Mnomain>)
endif ()
