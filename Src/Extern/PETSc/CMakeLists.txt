if (DIM EQUAL 1)
   message(FATAL_ERROR "PETSc interfaces are not supported for 1D builds")
endif ()

find_package(PETSc 2.13 REQUIRED)

target_compile_definitions( amrex
   PUBLIC
   $<BUILD_INTERFACE:AMREX_USE_PETSC>)

target_include_directories( amrex
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( amrex
   PRIVATE
   AMReX_PETSc.cpp
   AMReX_PETSc.H
   AMReX_PETSc_fort_mod.F90
   )

target_link_libraries( amrex PUBLIC PETSC )
