add_amrex_define(AMREX_USE_OPENPMD NO_LEGACY IF AMReX_OPENPMD)

foreach(D IN LISTS AMReX_SPACEDIM)

  target_include_directories( amrex_${D}d
     PUBLIC
     $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

  target_sources( amrex_${D}d
     PRIVATE
     AMReX_PlotFileUtilOPENPMD.H
     AMReX_PlotFileUtilOPENPMD.cpp
     AMReX_PlotFileOPENPMD.cpp
     AMReX_PlotFileOPENPMD_PTL.cpp
     AMReX_ParticlesOPENPMD.H
     AMReX_PlotFileUtilOPENPMD_PTLImpl.H
 )
endforeach()
