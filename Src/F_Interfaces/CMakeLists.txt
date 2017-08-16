set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")

set ( CXXSRC
   AmrCore/AMReX_amrcore_fi.cpp AmrCore/AMReX_tagbox_fi.cpp AmrCore/AMReX_FAmrCore.cpp
   AmrCore/AMReX_fillpatch_fi.cpp AmrCore/AMReX_fluxregister_fi.cpp
   Base/AMReX_boxarray_fi.cpp   Base/AMReX_multifab_fi.cpp
   Base/AMReX_plotfile_fi.cpp   Base/AMReX_box_fi.cpp       Base/AMReX_FPhysBC.cpp
   Base/AMReX_multifabutil_fi.cpp  Base/AMReX_vismf_fi.cpp  Base/AMReX_distromap_fi.cpp
   Base/AMReX_geometry_fi.cpp   Base/AMReX_physbc_fi.cpp
   Octree/AMReX_octree_fi.cpp )

set ( F90SRC
   AmrCore/AMReX_amrcore_mod.F90       AmrCore/AMReX_fillpatch_mod.F90
   AmrCore/AMReX_interpolater_mod.F90  AmrCore/AMReX_amr_mod.F90
   AmrCore/AMReX_fluxregister_mod.F90  AmrCore/AMReX_tagbox_mod.F90
   Base/AMReX_base_mod.F90      Base/AMReX_distromap_mod.F90     Base/AMReX_init_mod.F90
   Base/AMReX_parallel_mod.F90  Base/AMReX_boxarray_mod.F90      Base/AMReX_fab_mod.F90
   Base/AMReX_multifab_mod.F90  Base/AMReX_physbc_mod.F90        Base/AMReX_box_mod.F90
   Base/AMReX_geometry_mod.F90  Base/AMReX_multifabutil_mod.F90  Base/AMReX_plotfile_mod.F90
   Octree/AMReX_octree_mod.F90 )

set ( ALLHEADERS AmrCore/AMReX_FAmrCore.H Base/AMReX_FPhysBC.H )


# Accumulate sources
set ( ALLSRC ${CXXSRC} ${F90SRC} )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC     ${CMAKE_CURRENT_LIST_DIR} )
prepend ( ALLHEADERS ${CMAKE_CURRENT_LIST_DIR} )

target_sources ( amrex PRIVATE ${ALLSRC} ${ALLHEADERS} )

# 
# Accumulate headers 
#
accumulate_headers (ALLHEADERS)   

   
