# -*- mode: cmake -*-------------------------------------------
cmake_minimum_required(VERSION 2.8)

if (ENABLE_MPI)
  set(CMAKE_C_COMPILER ${MPI_PREFIX}/bin/mpicc)
  set(CMAKE_CXX_COMPILER ${MPI_PREFIX}/bin/mpic++)
  set(CMAKE_Fortran_COMPILER ${MPI_PREFIX}/bin/mpif90)
endif()

enable_language(C)
enable_language(CXX)
enable_language(Fortran)

set(BL_MACHINE ${CMAKE_SYSTEM_NAME})

message(STATUS "   BL_SPACEDIM = ${BL_SPACEDIM} (1,2,3)")
message(STATUS "   BL_MACHINE = ${BL_MACHINE} (<ARCH>)")
message(STATUS "   BL_PRECISION = ${BL_PRECISION} (FLOAT, DOUBLE)")
message(STATUS "   ENABLE_MPI = ${ENABLE_MPI} (0,1)")
message(STATUS "   ENABLE_OpenMP = ${ENABLE_OpenMP} (0,1)")

if (ENABLE_OpenMP)
  find_package(OpenMP)
  list(APPEND CMAKE_CC_FLAGS "${OpenMP_C_FLAGS}")
  list(APPEND CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS}")
endif(ENABLE_OpenMP)

include(FortranCInterface)
include(${FortranCInterface_BINARY_DIR}/Output.cmake)

message(STATUS "BoxLib-specific compile settings:")
if(FortranCInterface_GLOBAL_SUFFIX STREQUAL ""  AND FortranCInterface_GLOBAL_CASE STREQUAL "UPPER")
    message(STATUS "   Fortran name mangling scheme to UPPERCASE (upper case, no append underscore)")
    set(BL_FORTLINK UPPERCASE)
elseif(FortranCInterface_GLOBAL_SUFFIX STREQUAL ""  AND FortranCInterface_GLOBAL_CASE STREQUAL "LOWER")
    message(STATUS "   Fortran name mangling scheme to LOWERCASE (lower case, no append underscore)")
    set(BL_FORTLINK LOWERCASE)
elseif(FortranCInterface_GLOBAL_SUFFIX STREQUAL "_" AND FortranCInterface_GLOBAL_CASE STREQUAL "LOWER")
    message(STATUS "   Fortran name mangling scheme to UNDERSCORE (lower case, append underscore)")
    set(BL_FORTLINK "UNDERSCORE")
else()
    message(AUTHOR_WARNING "Fortran to C mangling not backward compatible with older style BoxLib code") 
endif()

set(BL_DEFINES BL_NOLINEVALUES;BL_PARALLEL_IO;BL_SPACEDIM=${BL_SPACEDIM};BL_FORT_USE_${BL_FORTLINK};BL_${BL_MACHINE};BL_USE_${BL_PRECISION})

if (ENABLE_MPI)
  set(MPI_SUFFIX .MPI)
  list(APPEND BL_DEFINES BL_USE_MPI)
else()
  set(MPI_SUFFIX)
endif()

if (ENABLE_OpenMP)
  set(OMP_SUFFIX .OMP)
  find_package(OpenMP REQUIRED)
  list(APPEND BL_DEFINES BL_USE_OMP)
else()
  set(OMP_SUFFIX)
endif()

