set (CXXSRC "")
set (F90SRC "")
set (F77SRC "")
set (ALLHEADERS "")
set (ALLSRC "")

set ( CXXSRC AMReX_TracerParticles.cpp    AMReX_LoadBalanceKD.cpp)
   
set ( F90SRC AMReX_Particle_mod_${BL_SPACEDIM}d.F90    AMReX_KDTree_${BL_SPACEDIM}d.F90)

list ( APPEND F90SRC AMReX_OMPDepositionHelper_nd.F90 )

set ( F77SRC AMReX_Particles_${BL_SPACEDIM}D.F )

set ( ALLHEADERS 
   AMReX_NeighborParticles.H   AMReX_ParGDB.H    AMReX_ParticleContainerI.H
   AMReX_ParticleInit.H  AMReX_Particles.H   AMReX_NeighborParticlesI.H
   AMReX_ParIterI.H  AMReX_ParticleI.H    AMReX_Particles_F.H
   AMReX_TracerParticles.H    AMReX_LoadBalanceKD.H    AMReX_KDTree_F.H)

# Accumulate sources
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
