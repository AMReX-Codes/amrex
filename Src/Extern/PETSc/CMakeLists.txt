if (AMReX_SPACEDIM EQUAL 1)
   message(FATAL_ERROR "PETSc interfaces are not supported for 1D builds")
endif ()

add_amrex_define(AMREX_USE_PETSC NO_LEGACY)

foreach(D IN LISTS AMReX_SPACEDIM)
    target_include_directories(amrex_${D}d
       PUBLIC
       $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

    target_sources(amrex_${D}d
       PRIVATE
       AMReX_PETSc.cpp
       AMReX_PETSc.H
    )
endforeach()
