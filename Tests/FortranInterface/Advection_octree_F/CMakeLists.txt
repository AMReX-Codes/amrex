foreach(D IN LISTS AMReX_SPACEDIM)
    if (NOT D EQUAL 2)
       continue()
    endif ()

    set(_sources advect_2d_mod.F90 compute_flux_2d.f90 slope_2d.f90)
    list(TRANSFORM _sources PREPEND Src_2d/ )

    list(APPEND _sources amr_data_mod.F90 averagedown_mod.F90 bc_mod.F90
       compute_dt_mod.F90 evolve_mod.F90 fillpatch_mod.F90 fmain.F90
       initdata.F90 my_amr_mod.F90 plotfile_mod.F90 tagging_mod.F90)

    list(TRANSFORM _sources PREPEND Source/ )

    list(APPEND _sources
        Exec/SingleVortex/face_velocity_2d.F90
        Exec/SingleVortex/Prob.f90
    )

    set(_input_files inputs.rt)
    list(TRANSFORM _input_files PREPEND "Exec/SingleVortex/")

    setup_test(${D} _sources _input_files
        CMDLINE_PARAMS_AFTER_INPUTS max_step=5
        HAS_FORTRAN_MODULES)

    unset(_sources)
    unset(_input_files)
endforeach()
