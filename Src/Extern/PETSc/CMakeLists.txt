if (AMReX_SPACEDIM EQUAL 1)
   message(FATAL_ERROR "PETSc interfaces are not supported for 1D builds")
endif ()

add_amrex_define(AMREX_USE_PETSC NO_LEGACY)

target_include_directories( amrex
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( amrex
   PRIVATE
   AMReX_PETSc.cpp
   AMReX_PETSc.H
   )
