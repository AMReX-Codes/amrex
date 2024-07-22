add_amrex_define(AMREX_USE_PETSC NO_LEGACY NO_1D IF AMReX_PETSC)

foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
       message(STATUS "PETSc interfaces are not supported for 1D builds (skipping)")
       continue()
    endif ()

    target_include_directories(amrex_${D}d
       PUBLIC
       $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

    target_sources(amrex_${D}d
       PRIVATE
       AMReX_PETSc.cpp
       AMReX_PETSc.H
    )
endforeach()
