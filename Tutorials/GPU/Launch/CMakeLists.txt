#
# This test requires CUDA to be enabled
#
if (NOT ENABLE_CUDA)
   return ()
endif ()

set( SRC_DIR ${CMAKE_CURRENT_LIST_DIR} )
set( EXENAME  "gpu_launch.exe" )

add_executable( ${EXENAME} EXCLUDE_FROM_ALL "")

target_sources( ${EXENAME}
   PRIVATE
   ${SRC_DIR}/MyKernel.H
   ${SRC_DIR}/MyKernel_F.F90
   ${SRC_DIR}/MyKernel_F.H
   ${SRC_DIR}/main.cpp)

set_source_files_properties(${SRC_DIR}/main.cpp PROPERTIES LANGUAGE CUDA)

get_target_property(amrex_flags  amrex COMPILE_OPTIONS)

evaluate_genex(amrex_flags amrex_fortran_flags
   LANG   Fortran
   COMP   ${CMAKE_Fortran_COMPILER_ID}
   CONFIG ${CMAKE_BUILD_TYPE}
   STRING )

# Since we are forcing the use of fortran compiler to link
# we need to specify the flags to add at link phase since
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
   "${amrex_fortran_flags}"
   CUDA_STANDARD 11               # Adds -std=c++11
   CUDA_STANDARD_REQUIRED ON
   CUDA_SEPARABLE_COMPILATION ON  
   CUDA_RESOLVE_DEVICE_SYMBOLS OFF # This is crucial to avoid the device code linkage step
   )

target_link_libraries(${EXENAME} amrex)

add_tutorial(${EXENAME})   
