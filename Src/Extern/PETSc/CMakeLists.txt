if (DIM EQUAL 1)
   message(FATAL_ERROR "PETSc interfaces are not supported for 1D builds")
endif ()

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
   )
