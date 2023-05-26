foreach(D IN LISTS AMReX_SPACEDIM)
    target_include_directories(amrex_${D}d PUBLIC
       $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

    target_sources(amrex_${D}
       PRIVATE
       AMReX_DG.H
       AMReX_DG.cpp
       )

   if (AMReX_FORTRAN)
       target_sources(amrex
          PRIVATE
          AMReX_DGInterfaceModule.F90
          AMReX_DGInterfaceModule_fi.cpp
          )
    endif ()
endforeach()
