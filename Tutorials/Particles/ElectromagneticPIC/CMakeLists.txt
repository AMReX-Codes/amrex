#
# Contains 3 tutorials
#
#  * CUDA
#  * OpenMP
#  * OpenACC (NOT supported anymore)
#
if ( (NOT ENABLE_CUDA AND NOT ENABLE_OMP) OR (NOT DIM EQUAL 3) )
   return()
endif ()


set(_sources Constants.H EMParticleContainer.H EMParticleContainerInit.cpp Evolve.H
   IO.cpp IO.H main.cpp NodalFlags.cpp NodalFlags.H )

list(TRANSFORM _sources PREPEND Source/)


##########################################################################################
#
# CUDA tutorial  --- This does not work
#
##########################################################################################
# if (ENABLE_CUDA)
#    set(_cuda_sources EMParticleContainer.cpp  em_pic_K.H  Evolve.cpp)
#    list(TRANSFORM _cuda_sources PREPEND Exec/CUDA/)
#    list(APPEND _cuda_sources ${_sources})

#    set(_input_files Exec/CUDA/inputs)

#    setup_tutorial(_cuda_sources _input_files
#       BASE_NAME ElectromagneticPIC_CUDA
#       RUNTIME_SUBDIR CUDA)

#    unset(_cuda_sources)
# endif ()


##########################################################################################
#
# OpenMP tutorial
#
##########################################################################################
if (ENABLE_OMP)
   set(_omp_sources EMParticleContainer.cpp em_pic_3d.F90 em_pic_F.H Evolve.cpp)
   list(TRANSFORM _omp_sources PREPEND Exec/OpenMP/)
   list(APPEND _omp_sources ${_sources})

   set(_input_files Exec/OpenMP/inputs)

   setup_tutorial(_omp_sources _input_files
      BASE_NAME ElectromagneticPIC_OpenMP
      RUNTIME_SUBDIR OpenMP)

   unset(_omp_sources)
endif ()


unset(_sources)
unset(_input_files)
