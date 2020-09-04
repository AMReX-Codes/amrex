if ( NOT ENABLE_PARTICLES OR (DIM EQUAL 1) )
   return ()
endif ()

message(WARNING "\nTutorial FortranInterfaces/Advection_F fails at runtime if ENABLE_EB=ON")

set(_sources advect_${DIM}d_mod.F90 compute_flux_${DIM}d.f90 slope_${DIM}d.f90 )
list(TRANSFORM _sources PREPEND Src_${DIM}d/ )

list(APPEND _sources amr_data_mod.F90 averagedown_mod.F90 bc_mod.F90 compute_dt_mod.F90
   evolve_mod.F90 fillpatch_mod.F90 fmain.F90 initdata.F90 Make.package my_amr_mod.F90
   plotfile_mod.F90 tagging_mod.F90)

list(TRANSFORM _sources PREPEND Source/ )

list(APPEND _sources
   Exec/SingleVortex/face_velocity_${DIM}d.F90
   Exec/SingleVortex/Prob_${DIM}d.f90
   )

file( GLOB_RECURSE _input_files LIST_DIRECTORIES false
   ${CMAKE_CURRENT_LIST_DIR}/input* )

setup_tutorial(_sources _input_files HAS_FORTRAN_MODULES)

unset(_sources)
unset(_input_files)
