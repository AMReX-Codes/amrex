
#
# Compile definitions required by codes using Fortran interfaces
# These definitions are not present anywhere in AMReX source code
# so I think they are meant to be used by the app code only
#
add_amrex_define(AMREX_USE_F_INTERFACES)

#
# BASE subdir
#
target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/Base>)

target_sources( amrex PRIVATE
   Base/AMReX_init_mod.F90
   Base/AMReX_parallel_mod.F90
   Base/AMReX_base_mod.F90
   Base/AMReX_box_mod.F90
   Base/AMReX_boxarray_mod.F90
   Base/AMReX_distromap_mod.F90
   Base/AMReX_geometry_mod.F90
   Base/AMReX_multifab_mod.F90
   Base/AMReX_plotfile_mod.F90
   Base/AMReX_multifabutil_mod.F90
   Base/AMReX_physbc_mod.F90
   Base/AMReX_fab_mod.F90
   Base/AMReX_init_fi.cpp
   Base/AMReX_box_fi.cpp
   Base/AMReX_boxarray_fi.cpp
   Base/AMReX_distromap_fi.cpp
   Base/AMReX_geometry_fi.cpp
   Base/AMReX_multifab_fi.cpp
   Base/AMReX_vismf_fi.cpp
   Base/AMReX_plotfile_fi.cpp
   Base/AMReX_multifabutil_fi.cpp
   Base/AMReX_physbc_fi.cpp
   Base/AMReX_FPhysBC.H
   Base/AMReX_FPhysBC.cpp
   )

if (AMReX_MPI)
  target_sources( amrex PRIVATE
    Base/AMReX_fi_mpi_mod.F90
    Base/AMReX_mpi_reduce_int.F90
    Base/AMReX_mpi_reduce_real.F90
    )
endif ()

#
# AMRCORE subdir
#
target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/AmrCore>)

target_sources(amrex PRIVATE
   AmrCore/AMReX_amrcore_fi.cpp
   AmrCore/AMReX_amrcore_mod.F90
   AmrCore/AMReX_amr_mod.F90
   AmrCore/AMReX_FAmrCore.cpp
   AmrCore/AMReX_FAmrCore.H
   AmrCore/AMReX_fillpatch_fi.cpp
   AmrCore/AMReX_fillpatch_mod.F90
   AmrCore/AMReX_FlashFluxRegister.cpp
   AmrCore/AMReX_flash_fluxregister_fi.cpp
   AmrCore/AMReX_FlashFluxRegister.H
   AmrCore/AMReX_flash_fluxregister_mod.F90
   AmrCore/AMReX_fluxregister_fi.cpp
   AmrCore/AMReX_fluxregister_mod.F90
   AmrCore/AMReX_interpolater_mod.F90
   AmrCore/AMReX_tagbox_fi.cpp
   AmrCore/AMReX_tagbox_mod.F90
   )


#
# OCTREE subdir
#
target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/Octree>)

target_sources( amrex PRIVATE Octree/AMReX_octree_mod.F90  Octree/AMReX_octree_fi.cpp )


#
# LINEARSOLVERS subdir
#
if (AMReX_LINEAR_SOLVERS)
   target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/LinearSolvers>)

   target_sources( amrex PRIVATE
      LinearSolvers/AMReX_abeclaplacian_fi.cpp
      LinearSolvers/AMReX_abeclaplacian_mod.F90
      LinearSolvers/AMReX_linear_solver_mod.F90
      LinearSolvers/AMReX_linop_mod.F90
      LinearSolvers/AMReX_multigrid_mod.F90
      LinearSolvers/AMReX_multigrid_fi.cpp
      LinearSolvers/AMReX_poisson_mod.F90
      LinearSolvers/AMReX_linop_fi.cpp
      LinearSolvers/AMReX_poisson_fi.cpp
      )
endif ()

#
# PARTICLE subdir
#
if ( AMReX_PARTICLES )
   target_include_directories(amrex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/Particle>)
   target_sources( amrex PRIVATE
      Particle/AMReX_particlecontainer_fi.cpp Particle/AMReX_particlecontainer_mod.F90 )
endif ()
