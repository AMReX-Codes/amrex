#
# This test requires CUDA to be enabled
#
if (NOT ENABLE_CUDA)
   return ()
endif ()

set( SRC_DIR ${CMAKE_CURRENT_LIST_DIR} )
set( EXENAME  "gpu_launch.exe" )

print(EXENAME)


#
# If ENABLE_CUDA_FORTRAN is ON, this test cannot be compiled yes
# because of lack of support of CUDA C++/Fortran in Cmake 3.11
#
if (ENABLE_CUDA_FORTRAN)
   message(WARNING "Compilation of ${EXENAME} with ENABLE_CUDA_FORTRAN=ON "
      "is not yet supported. Skipping test")
   return()
endif ()


add_executable( ${EXENAME} EXCLUDE_FROM_ALL "")

target_sources( ${EXENAME}
   PRIVATE
   ${SRC_DIR}/MyKernel.H
   ${SRC_DIR}/MyKernel_F.F90
   ${SRC_DIR}/MyKernel_F.H
   ${SRC_DIR}/main.cpp)

set_source_files_properties(${SRC_DIR}/main.cpp PROPERTIES LANGUAGE CUDA)


#if (ENABLE_FORTRAN_CUDA)
   string(REPLACE ";" " " tmp "${NVCC_ARCH_FLAGS}")
   set(CMAKE_CUDA_FLAGS "${tmp}")

   # Since we are forcing the use of fortran compiler to link
   # we need to specify the flags to add since
   # it won't propagate amrex ones
   set_target_properties( ${EXENAME}
      PROPERTIES
      Fortran_MODULE_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR}/mod_files
      RUNTIME_OUTPUT_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR}/Launch
      LINKER_LANGUAGE
      Fortran
      LINK_FLAGS
      "-Mcuda=cc70,ptxinfo,fastmath,charstring -Mnomain -noacc CUDA_HOME=/sw/summit/cuda/9.2.148"
      CUDA_STANDARD 11               # Adds -std=c++11
      CUDA_STANDARD_REQUIRED ON
      CUDA_SEPARABLE_COMPILATION ON
      CUDA_RESOLVE_DEVICE_SYMBOLS OFF
      # LINK_FLAGS
      # "-Mnomain -Mcuda=cc70,ptxinfo,fastmath,charstring"
      )

# else ()
   
#    set_target_properties( ${EXENAME}
#       PROPERTIES
#       Fortran_MODULE_DIRECTORY
#       ${CMAKE_CURRENT_BINARY_DIR}/mod_files
#       RUNTIME_OUTPUT_DIRECTORY
#       ${CMAKE_CURRENT_BINARY_DIR}/Launch
#       LINKER_LANGUAGE
#       CUDA
#       LINK_FLAGS
#       -std=c++11
#       CXX_STANDARD 11
#       CXX_STANDARD_REQUIRED ON
#       CUDA_STANDARD 11               # Adds -std=c++11
#       CUDA_STANDARD_REQUIRED ON
#       CUDA_SEPARABLE_COMPILATION ON
#       CUDA_RESOLVE_DEVICE_SYMBOLS OFF
#       )   
   
# endif ()

# set_target_properties(${EXENAME} PROPERTIES CUDA_SEPARABLE_COMPILATION OFF)

target_link_libraries(${EXENAME} amrex cuda)

add_tutorial(${EXENAME})   
