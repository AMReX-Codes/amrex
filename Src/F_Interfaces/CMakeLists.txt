#---------------------------------------------#
#  BASE subdir
#---------------------------------------------#
set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")

list ( APPEND F90SRC      AMReX_init_mod.F90 AMReX_parallel_mod.F90 AMReX_base_mod.F90 )
list ( APPEND F90SRC      AMReX_box_mod.F90 AMReX_boxarray_mod.F90 AMReX_distromap_mod.F90 )
list ( APPEND F90SRC      AMReX_geometry_mod.F90 AMReX_multifab_mod.F90 AMReX_plotfile_mod.F90 )
list ( APPEND F90SRC      AMReX_multifabutil_mod.F90 AMReX_physbc_mod.F90 AMReX_fab_mod.F90 )

list ( APPEND CXXSRC      AMReX_init_fi.cpp AMReX_box_fi.cpp AMReX_boxarray_fi.cpp AMReX_distromap_fi.cpp )
list ( APPEND CXXSRC      AMReX_geometry_fi.cpp AMReX_multifab_fi.cpp AMReX_vismf_fi.cpp AMReX_plotfile_fi.cpp )
list ( APPEND CXXSRC      AMReX_multifabutil_fi.cpp AMReX_physbc_fi.cpp )

list ( APPEND ALLHEADERS  AMReX_FPhysBC.H )
list ( APPEND CXXSRC      AMReX_FPhysBC.cpp )

#
# Accumulate sources
#
set ( ALLSRC ${CXXSRC} ${F90SRC} )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC     ${CMAKE_CURRENT_LIST_DIR}/Base )
prepend ( ALLHEADERS ${CMAKE_CURRENT_LIST_DIR}/Base )

target_sources ( amrex PRIVATE ${ALLSRC} ${ALLHEADERS} )

# 
# Accumulate headers 
#
accumulate_headers (ALLHEADERS)   

#
# Add to include dirs
# 
accumulate_include_dirs (${CMAKE_CURRENT_LIST_DIR}/Base)   


#---------------------------------------------#
# AMRCORE subdir
#---------------------------------------------#
set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")

list ( APPEND F90SRC      AMReX_amr_mod.F90 AMReX_amrcore_mod.F90 AMReX_tagbox_mod.F90 )
list ( APPEND F90SRC      AMReX_fillpatch_mod.F90 AMReX_fluxregister_mod.F90 AMReX_interpolater_mod.F90 )

list ( APPEND CXXSRC      AMReX_amrcore_fi.cpp AMReX_tagbox_fi.cpp AMReX_fillpatch_fi.cpp )
list ( APPEND CXXSRC      AMReX_fluxregister_fi.cpp AMReX_FAmrCore.cpp )

list ( APPEND ALLHEADERS  AMReX_FAmrCore.H )


#
# Accumulate sources
#
set ( ALLSRC ${CXXSRC} ${F90SRC} )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC     ${CMAKE_CURRENT_LIST_DIR}/AmrCore )
prepend ( ALLHEADERS ${CMAKE_CURRENT_LIST_DIR}/AmrCore )

target_sources ( amrex PRIVATE ${ALLSRC} ${ALLHEADERS} )

# 
# Accumulate headers 
#
accumulate_headers (ALLHEADERS)   

#
# Add to include dirs
# 
accumulate_include_dirs (${CMAKE_CURRENT_LIST_DIR}/AmrCore)   


#---------------------------------------------#
# OCTREE subdir
#---------------------------------------------#
set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")

list ( APPEND F90SRC    AMReX_octree_mod.F90 )

list ( APPEND CXXSRC    AMReX_octree_fi.cpp )


#
# Accumulate sources
#
set ( ALLSRC ${CXXSRC} ${F90SRC} )

#
# Add absolute path (needed by target_sources)
# 
prepend ( ALLSRC     ${CMAKE_CURRENT_LIST_DIR}/Octree )
prepend ( ALLHEADERS ${CMAKE_CURRENT_LIST_DIR}/Octree )

target_sources ( amrex PRIVATE ${ALLSRC} ${ALLHEADERS} )

# 
# Accumulate headers 
#
accumulate_headers (ALLHEADERS)   

#
# Add to include dirs
# 
accumulate_include_dirs (${CMAKE_CURRENT_LIST_DIR}/Octree)   
