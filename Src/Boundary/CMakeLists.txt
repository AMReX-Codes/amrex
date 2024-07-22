foreach(D IN LISTS AMReX_SPACEDIM)
    target_include_directories(amrex_${D}d PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

    target_sources(amrex_${D}d
       PRIVATE
       AMReX_Mask.cpp
       AMReX_MultiMask.cpp
       AMReX_FabSet.H
       AMReX_BndryRegister.H
       AMReX_Mask.H
       AMReX_MultiMask.H
       AMReX_BndryData.H
       AMReX_BoundCond.H
       AMReX_InterpBndryData.H
       AMReX_LO_BCTYPES.H
       AMReX_LO_BCTYPES.cpp
       AMReX_InterpBndryData_K.H
       AMReX_InterpBndryData_${D}D_K.H
       AMReX_LOUtil_K.H
       AMReX_YAFluxRegister.H
       AMReX_YAFluxRegister_K.H
       AMReX_YAFluxRegister_${D}D_K.H
       AMReX_BoundaryFwd.H
       )

    if (NOT (D EQUAL 1))
       target_sources(amrex_${D}d
          PRIVATE
          AMReX_EdgeFluxRegister.H
          AMReX_EdgeFluxRegister.cpp
          )
    endif ()

    if (AMReX_FORTRAN)
       target_sources(amrex_${D}d
          PRIVATE
          AMReX_LO_UTIL.F90
          AMReX_lo_bctypes_mod.F90
          )
    endif ()
endforeach()
