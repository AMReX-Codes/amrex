if (NOT AMReX_AMRLEVEL)
    return()
endif ()

foreach(D IN LISTS AMReX_SPACEDIM)
    set(_sources main.cpp WMLES_advance.cpp WMLES_advance_box.cpp WMLES_advance_box_eb.cpp WMLES_bcfill.cpp
       WMLESBld.cpp WMLES.cpp WMLES.H WMLES_derive.cpp WMLES_derive.H WMLES_index_macros.H WMLES_init_eb2.cpp WMLES_io.cpp
       WMLES_K.H WMLES_parm.cpp WMLES_parm.H WMLES_setup.cpp WMLES_tagging.H
       WMLES_sgs_K.H WMLES_wallmodel_K.H)

    set(_diffusion_sources WMLES_diffusion_K.H WMLES_diffusion_eb_K.H )
    list(TRANSFORM _diffusion_sources PREPEND diffusion/)

    set(_hydro_sources WMLES_divop_K.H WMLES_hydro_eb_K.H WMLES_hydro_K.H )
    list(TRANSFORM _hydro_sources PREPEND hydro/)

    list(APPEND _sources ${_diffusion_sources} ${_hydro_sources})
    list(TRANSFORM _sources PREPEND Source/)

    unset(_diffusion_sources)
    unset(_hydro_sources)

    ##########################################################################
    #
    # Cylinder test — crossflow over a cylinder with WMLES
    #
    ##########################################################################
    set(_cylinder_sources Exec/Cylinder/wmles_prob.H  Exec/Cylinder/wmles_prob.cpp
                          Exec/Cylinder/wmles_prob_parm.H  Exec/Cylinder/wmles_prob_parm.cpp
                          Exec/Cylinder/WMLES_bcfill.cpp ${_sources})
    set(_input_files Exec/Cylinder/inputs Exec/Cylinder/inputs-ci)

    setup_test(${D} _cylinder_sources _input_files
       BASE_NAME WMLES_Cylinder
       RUNTIME_SUBDIR Cylinder)

    unset(_cylinder_sources)

    ##########################################################################
    #
    # STL_Box test — 3D flow over an STL-imported box (3D only)
    #
    ##########################################################################
    if (D EQUAL 3)
        set(_stlbox_sources Exec/STL_Box/wmles_prob.H  Exec/STL_Box/wmles_prob.cpp
                            Exec/STL_Box/wmles_prob_parm.H  Exec/STL_Box/wmles_prob_parm.cpp
                            Exec/STL_Box/WMLES_bcfill.cpp ${_sources})
        set(_stlbox_inputs Exec/STL_Box/inputs Exec/STL_Box/inputs-ci)

        setup_test(${D} _stlbox_sources _stlbox_inputs
           BASE_NAME WMLES_STL_Box
           RUNTIME_SUBDIR STL_Box)

        unset(_stlbox_sources)
        unset(_stlbox_inputs)
    endif()

    # Clean-up
    unset(_sources)
    unset(_input_files)
endforeach()
