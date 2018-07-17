#---------------------------------------------#
# BASE subdir
#---------------------------------------------#
add_sources ( Base/AMReX_init_mod.F90         Base/AMReX_parallel_mod.F90 Base/AMReX_base_mod.F90 )
add_sources ( Base/AMReX_box_mod.F90          Base/AMReX_boxarray_mod.F90 Base/AMReX_distromap_mod.F90 )
add_sources ( Base/AMReX_geometry_mod.F90     Base/AMReX_multifab_mod.F90 Base/AMReX_plotfile_mod.F90 )
add_sources ( Base/AMReX_multifabutil_mod.F90 Base/AMReX_physbc_mod.F90   Base/AMReX_fab_mod.F90 )

add_sources ( Base/AMReX_init_fi.cpp         Base/AMReX_box_fi.cpp      Base/AMReX_boxarray_fi.cpp Base/AMReX_distromap_fi.cpp )
add_sources ( Base/AMReX_geometry_fi.cpp     Base/AMReX_multifab_fi.cpp Base/AMReX_vismf_fi.cpp    Base/AMReX_plotfile_fi.cpp )
add_sources ( Base/AMReX_multifabutil_fi.cpp Base/AMReX_physbc_fi.cpp )

add_sources ( Base/AMReX_FPhysBC.H )
add_sources ( Base/AMReX_FPhysBC.cpp )


#---------------------------------------------#
# AMRCORE subdir
#---------------------------------------------#
add_sources ( AmrCore/AMReX_amr_mod.F90       AmrCore/AMReX_amrcore_mod.F90      AmrCore/AMReX_tagbox_mod.F90 )
add_sources ( AmrCore/AMReX_fillpatch_mod.F90 AmrCore/AMReX_fluxregister_mod.F90 AmrCore/AMReX_interpolater_mod.F90 )

add_sources ( AmrCore/AMReX_amrcore_fi.cpp      AmrCore/AMReX_tagbox_fi.cpp AmrCore/AMReX_fillpatch_fi.cpp )
add_sources ( AmrCore/AMReX_fluxregister_fi.cpp AmrCore/AMReX_FAmrCore.cpp )

add_sources ( AmrCore/AMReX_FAmrCore.H )



#---------------------------------------------#
# OCTREE subdir
#---------------------------------------------#
add_sources ( Octree/AMReX_octree_mod.F90 )
add_sources ( Octree/AMReX_octree_fi.cpp )


#---------------------------------------------#
# LINEARSOLVERS subdir
#---------------------------------------------#
if (ENABLE_LINEAR_SOLVERS)
   add_sources ( LinearSolvers/AMReX_abeclaplacian_fi.cpp   LinearSolvers/AMReX_abeclaplacian_mod.F90 )
   add_sources ( LinearSolvers/AMReX_linear_solver_mod.F90  LinearSolvers/AMReX_linop_mod.F90 )
   add_sources ( LinearSolvers/AMReX_multigrid_mod.F90      LinearSolvers/AMReX_multigrid_fi.cpp )
   add_sources ( LinearSolvers/AMReX_poisson_mod.F90        LinearSolvers/AMReX_linop_fi.cpp )
   add_sources ( LinearSolvers/AMReX_poisson_fi.cpp )    
endif ()

#---------------------------------------------#
# PARTICLE subdir
#---------------------------------------------#
if ( ENABLE_PARTICLES )
   add_sources ( Particle/AMReX_particlecontainer_fi.cpp Particle/AMReX_particlecontainer_mod.F90 )
endif ()
