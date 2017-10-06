set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")


list ( APPEND CXXSRC      AMReX_TracerParticles.cpp AMReX_LoadBalanceKD.cpp)

list ( APPEND ALLHEADERS  AMReX_Particles.H AMReX_ParGDB.H AMReX_TracerParticles.H )
list ( APPEND ALLHEADERS  AMReX_NeighborParticles.H AMReX_NeighborParticlesI.H )
list ( APPEND ALLHEADERS  AMReX_ParticleI.H AMReX_ParticleInit.H AMReX_ParticleContainerI.H )
list ( APPEND ALLHEADERS  AMReX_LoadBalanceKD.H AMReX_KDTree_F.H )
list ( APPEND ALLHEADERS  AMReX_ParIterI.H AMReX_Particles_F.H )

list ( APPEND F77SRC      AMReX_Particles_${DIM}D.F )
list ( APPEND F90SRC      AMReX_Particle_mod_${DIM}d.F90 AMReX_KDTree_${DIM}d.F90)
list ( APPEND F90SRC      AMReX_OMPDepositionHelper_nd.F90 )

#
# Accumulate sources
# 
set ( ALLSRC ${CXXSRC} ${F90SRC} ${F77SRC} )

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


#
# Add to include dirs
# 
accumulate_include_dirs (${CMAKE_CURRENT_LIST_DIR})
