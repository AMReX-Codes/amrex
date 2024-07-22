foreach(D IN LISTS AMReX_SPACEDIM)
    target_include_directories(amrex_${D}d
       PUBLIC
       $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

    target_sources(amrex_${D}d
       PRIVATE
       AMReX_NVector_MultiFab.cpp
       AMReX_NVector_MultiFab.H
       AMReX_Sundials_Core.cpp
       AMReX_Sundials_Core.H
       AMReX_Sundials.H
       AMReX_SundialsIntegrator.H
       AMReX_SUNMemory.cpp
       AMReX_SUNMemory.H
    )
endforeach()
