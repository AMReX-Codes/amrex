target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(amrex
   PRIVATE
   AMReX_thornado.H
   AMReX_thornado.cpp
   )

if (AMReX_FORTRAN)
   target_sources(amrex
      PRIVATE
      AMReX_ThornadoInterfaceModule.F90
      AMReX_ThornadoInterfaceModule_fi.cpp
      )
endif ()
