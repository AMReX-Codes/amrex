target_include_directories( amrex
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( amrex
   PRIVATE
   AMReX_PlotFileUtilHDF5.H
   AMReX_PlotFileUtilHDF5.cpp
   AMReX_ParticleUtilHDF5.H
   AMReX_ParticleHDF5.H
   AMReX_ParticlesHDF5.H
   AMReX_WriteBinaryParticleDataHDF5.H
   )
