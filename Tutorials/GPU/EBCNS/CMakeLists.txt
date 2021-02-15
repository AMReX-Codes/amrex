if (NOT AMReX_SPACEDIM EQUAL 3 OR NOT AMReX_EB)
   return()
endif ()

set(_sources main.cpp CNS_advance.cpp CNS_bcfill.cpp CNSBld.cpp CNS.cpp CNS_derive.cpp CNS_derive.H CNS.H
   CNS_index_macros.H CNS_init_eb2.cpp CNS_io.cpp CNS_K.H CNS_parm.cpp CNS_parm.H CNS_setup.cpp CNS_tagging.H
   hydro/CNS_hydro_K.H)
list(TRANSFORM _sources PREPEND Source/)


set(_sources Exec/Sod/cns_prob.cpp Exec/Sod/cns_prob.H Exec/Sod/cns_prob_parm.H ${_sources})
set(_input_files Exec/Sod/inputs)

setup_tutorial(_sources _input_files
   BASE_NAME GPU_EBCNS_Sod
   RUNTIME_SUBDIR Sod
   EXTRA_DEFINITIONS AMREX_CRSEGRNDOMP)

unset(_sources)
unset(_input_files)
