# #
# # This test requires CUDA and Particles to be enabled
# #
# if (NOT AMReX_CUDA OR NOT AMReX_PARTICLES)
#    return ()
# endif ()
set(_sources CheckPair.H Constants.H main.cpp MDParticleContainer.cpp MDParticleContainer.H)
set(_input_files inputs)

setup_tutorial(_sources _input_files)

unset(_sources)
unset(_input_files)
