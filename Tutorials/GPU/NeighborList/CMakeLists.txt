#
# This test requires CUDA and Particles to be enabled
#
if (NOT ENABLE_CUDA OR NOT ENABLE_PARTICLES)
   return ()
endif ()

set( SRC_DIR ${CMAKE_CURRENT_LIST_DIR} )
set( EXENAME  "gpu_nblist.exe" )

add_executable( ${EXENAME} EXCLUDE_FROM_ALL "")

target_sources( ${EXENAME}
   PRIVATE
   ${SRC_DIR}/Constants.H
   ${SRC_DIR}/MDParticleContainer.H
   ${SRC_DIR}/MDParticleContainer.cpp
   ${SRC_DIR}/CheckPair.H
   ${SRC_DIR}/main.cpp)

target_include_directories(${EXENAME} PRIVATE ${CMAKE_CURRENT_LIST_DIR} )
set_source_files_properties(${SRC_DIR}/main.cpp ${SRC_DIR}/MDParticleContainer.cpp PROPERTIES LANGUAGE CUDA)

# Since we are forcing the use of fortran compiler to link
# we need to specify the flags to add at link phase since
# it won't propagate amrex ones
set_target_properties( ${EXENAME}
   PROPERTIES
   RUNTIME_OUTPUT_DIRECTORY
   ${CMAKE_CURRENT_BINARY_DIR}/Launch
   # LINK_FLAGS
   # "-Wl,--warn-unresolved-symbols"
   CXX_STANDARD 14               # Adds -std=c++11
   CXX_STANDARD_REQUIRED ON
   CUDA_STANDARD 14               # Adds -std=c++11
   CUDA_STANDARD_REQUIRED ON
   CUDA_SEPARABLE_COMPILATION ON  # This add -dc flag
   )

target_link_libraries(${EXENAME} amrex)

add_tutorial(${EXENAME})   
