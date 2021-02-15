if ( (AMReX_SPACEDIM EQUAL 1) OR NOT CMAKE_Fortran_COMPILER_LOADED )
   return()
endif ()


#
# This test works for both 2D and 3D builds
# There are two different tutorials here:
#
# 1) Single Vortex
# 2) Uniform Velocity
#


#
# Sources common to both cases
#
set(_sources Adv_F.H  AmrLevelAdv.cpp  AmrLevelAdv.H  LevelBldAdv.cpp  main.cpp)

list(APPEND _sources Src_${AMReX_SPACEDIM}d/Adv_${AMReX_SPACEDIM}d.f90 Src_${AMReX_SPACEDIM}d/compute_flux_${AMReX_SPACEDIM}d.f90
   Src_${AMReX_SPACEDIM}d/slope_${AMReX_SPACEDIM}d.f90 )

list(APPEND _sources Src_nd/Adv_nd.f90  Src_nd/Tagging_nd.f90  Src_nd/tagging_params.f90)

list(TRANSFORM _sources PREPEND Source/)


###############################################################################
#
# Single Vortex tutorial ------------------------------------------------------
#
###############################################################################
set(_sv_exe_dir Exec/SingleVortex/)

set(_sv_sources face_velocity_${AMReX_SPACEDIM}d.f90  Prob.f90)
list(TRANSFORM _sv_sources PREPEND ${_sv_exe_dir})
list(APPEND _sv_sources ${_sources})

set(_input_files inputs inputs.tracers probin)
list(TRANSFORM _input_files PREPEND ${_sv_exe_dir})

setup_tutorial(_sv_sources _input_files
   HAS_FORTRAN_MODULES
   BASE_NAME Advection_AmrLevel_SV
   RUNTIME_SUBDIR SingleVortex)

unset(_sv_sources)
unset(_sv_exe_dir)


###############################################################################
#
# Uniform Velocity tutorial ---------------------------------------------------
#
###############################################################################
set(_uv_exe_dir Exec/UniformVelocity/)

set(_uv_sources face_velocity_${AMReX_SPACEDIM}d.f90  Prob.f90 probdata.f90)
list(TRANSFORM _uv_sources PREPEND ${_uv_exe_dir})
list(APPEND _uv_sources ${_sources})

set(_input_files inputs inputs.regt probin)
list(TRANSFORM _input_files PREPEND ${_uv_exe_dir})

setup_tutorial(_uv_sources _input_files
   HAS_FORTRAN_MODULES
   BASE_NAME Advection_AmrLevel_UV
   RUNTIME_SUBDIR UniformVelocity)

unset(_uv_sources)
unset(_uv_exe_dir)


# Final clean up
unset(_sources)
unset(_input_files)
