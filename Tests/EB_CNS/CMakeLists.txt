
set(_sources main.cpp CNS_advance.cpp CNS_advance_box.cpp CNS_advance_box_eb.cpp CNS_bcfill.cpp
   CNSBld.cpp CNS.cpp CNS.H CNS_derive.cpp CNS_derive.H CNS_index_macros.H CNS_init_eb2.cpp CNS_io.cpp
   CNS_K.H CNS_parm.cpp CNS_parm.H CNS_setup.cpp CNS_tagging.H )

set(_diffusion_sources CNS_diffusion_K.H CNS_diffusion_eb_K.H )
list(TRANSFORM _diffusion_sources PREPEND diffusion/)

set(_hydro_sources CNS_divop_K.H CNS_flux_redistribute.cpp CNS_hydro_eb_K.H CNS_hydro_K.H )
list(TRANSFORM _hydro_sources PREPEND hydro/)

list(APPEND _sources ${_diffusion_sources} ${_hydro_sources})
list(TRANSFORM _sources PREPEND Source/)

unset(_diffusion_sources)
unset(_hydro_sources)

##########################################################################################
#
# Combustor tutorial
#
##########################################################################################
set(_combustor_sources Exec/Combustor/cns_prob.H  Exec/Combustor/cns_prob.cpp
                       Exec/Combustor/cns_prob_parm.H  Exec/Combustor/cns_prob_parm.cpp
                       Exec/Combustor/CNS_bcfill.cpp ${_sources})
set(_input_files Exec/Combustor/inputs Exec/Combustor/inputs-ci)

setup_test(_combustor_sources _input_files
   BASE_NAME CNS_Combustor
   RUNTIME_SUBDIR Combustor)

unset(_combustor_sources)

#
# Clean-up
#
unset(_sources)
unset(_input_files)
