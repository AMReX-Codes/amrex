#
# This directory contains 4 tutorials
#
#  * Combustor
#  * Pulse
#  * ShockRef
#  * Sod
#

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
set(_input_files Exec/Combustor/inputs Exec/Combustor/inputs.regt)

setup_test(_combustor_sources _input_files
   BASE_NAME CNS_Combustor
   RUNTIME_SUBDIR Combustor)

unset(_combustor_sources)



##########################################################################################
#
# Pulse tutorial
#
##########################################################################################
set(_pulse_sources Exec/Pulse/cns_prob.H Exec/Pulse/cns_prob.cpp
                   Exec/Pulse/cns_prob_parm.H ${_sources})
set(_input_files Exec/Pulse/inputs Exec/Pulse/inputs.regt)

setup_test(_pulse_sources _input_files
   BASE_NAME CNS_Pulse
   RUNTIME_SUBDIR Pulse)

unset(_pulse_sources)



##########################################################################################
#
# ShockRef tutorial
#
##########################################################################################
set(_pulse_sources Exec/ShockRef/cns_prob.H Exec/ShockRef/cns_prob.cpp
                   Exec/ShockRef/cns_prob_parm.H ${_sources})
set(_input_files Exec/ShockRef/inputs Exec/ShockRef/inputs.regt Exec/ShockRef/inputs.amr)

setup_test(_shockref_sources _input_files
   BASE_NAME CNS_ShockRef
   RUNTIME_SUBDIR ShockRef)

unset(_shockref_sources)



##########################################################################################
#
# Sod tutorial
#
##########################################################################################
set(_pulse_sources Exec/Sod/cns_prob.H Exec/Sod/cns_prob.cpp
                   Exec/Sod/cns_prob_parm.H ${_sources})
set(_input_files Exec/Sod/inputs)

setup_test(_sod_sources _input_files
   BASE_NAME CNS_Sod
   RUNTIME_SUBDIR Sod)

unset(_sod_sources)


#
# Clean-up
#
unset(_sources)
unset(_input_files)
