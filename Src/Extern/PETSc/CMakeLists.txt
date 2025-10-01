add_amrex_define(AMREX_USE_PETSC NO_LEGACY IF AMReX_PETSC)

foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
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
