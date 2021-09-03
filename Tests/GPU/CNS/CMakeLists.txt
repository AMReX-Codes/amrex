if (NOT AMReX_SPACEDIM EQUAL 3)
   return()
endif ()

set(_sources main.cpp CNS_advance.cpp CNS_bcfill.cpp CNSBld.cpp CNS.cpp CNS_derive.cpp CNS_derive.H CNS.H
   CNS_index_macros.H CNS_io.cpp CNS_K.H CNS_parm.cpp CNS_parm.H CNS_setup.cpp CNS_tagging.H
   hydro/CNS_hydro_K.H diffusion/CNS_diffusion_K.H)
list(TRANSFORM _sources PREPEND Source/)

##########################################################################################
#
# Sod tutorial
#
##########################################################################################
set(_sod_sources Exec/Sod/cns_prob.cpp Exec/Sod/cns_prob.H Exec/Sod/cns_prob_parm.H ${_sources})
set(_input_files Exec/Sod/inputs-ci)

setup_test(_sod_sources _input_files
   BASE_NAME GPU_CNS_Sod
   RUNTIME_SUBDIR Sod)

unset(_sod_sources)


#
# Clean-up
#
unset(_sources)
unset(_input_files)
