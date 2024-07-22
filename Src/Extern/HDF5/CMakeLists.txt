foreach(D IN LISTS AMReX_SPACEDIM)
    target_include_directories(amrex_${D}d
       PUBLIC
       $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

    target_sources(amrex_${D}d
       PRIVATE
       AMReX_PlotFileUtilHDF5.H
       AMReX_PlotFileUtilHDF5.cpp
       AMReX_ParticleUtilHDF5.H
       AMReX_ParticleHDF5.H
       AMReX_ParticlesHDF5.H
       AMReX_WriteBinaryParticleDataHDF5.H
    )
endforeach()
