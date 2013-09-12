# -*- mode: cmake -*-------------------------------------------
cmake_minimum_required(VERSION 2.8.8)
enable_language(CXX)
enable_language(C)
enable_language(Fortran)

# This is where the CCSE lib was installed
set(CCSE_DIR /home/marc/src/CCSE/JUNK)
set(HYPRE_DIR /home/marc/src/CASC/hyper)

# 1, 2, 3
set(BL_SPACEDIM 3)

# FLOAT DOUBLE
set(BL_PRECISION DOUBLE)

# 1 or 0 for both
set(ENABLE_MPI 1)
set(ENABLE_OpenMP 0)

# None, Debug, Release, RelWithDebInfo, MinSizeRel
set(CMAKE_BUILD_TYPE "Release")
set(EXTRA_DEFS_FOR_CCSE "CG_USE_OLD_CONVERGENCE_CRITERIA")


set(CMAKE_MODULE_PATH ${CCSE_DIR}/cmake)
include(CCSEOptions)
find_package(CCSE REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CCSE_INCLUDE_DIR})
link_directories(${CCSE_LIBRARY_DIR})

#
# Define a project name
# After this command the following varaibles are defined
#   MGC_SOURCE_DIR
#   MGC_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MGC)

set(CXX_sources  main.cpp writePlotFile.cpp)
set(CXX_includes writePlotFile.H)
set(F90_sources  )
set(F77_sources  )
set(FPP_includes COEF_F.H RHS_F.H)
set(FPP_sources  COEF_${BL_SPACEDIM}D.F RHS_${BL_SPACEDIM}D.F)

include(PreprocessBoxLibFortran)
preprocess_boxlib_fortran(FPP_out ${FPP_sources})

set(local_includes ${CXX_includes} ${FPP_includes})
set(local_sources  ${CXX_sources} ${F90_sources} ${F77_sources} ${FPP_out})

add_executable(mgc_tutorial ${local_sources} ${local_includes})
target_link_libraries(mgc_tutorial ${CCSE_LIBRARIES} ${MPI_LIBRARIES})

# Copy test directory files if an out of source build
if (NOT (${MGC_SOURCE_DIR} EQUAL ${MGC_BINARY_DIR}) )
  execute_process(COMMAND ${CMAKE_COMMAND} -E 
                  copy ${MGC_SOURCE_DIR}/inputs.3d ${MGC_BINARY_DIR}/inputs.3d) 
endif()
