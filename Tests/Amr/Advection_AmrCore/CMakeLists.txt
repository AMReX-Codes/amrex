foreach(D IN LISTS AMReX_SPACEDIM)
    if (D EQUAL 1)
       continue()
    endif ()

    # List of source files
    set(_sources AdvancePhiAllLevels.cpp  AdvancePhiAtLevel.cpp  AmrCoreAdv.cpp  AmrCoreAdv.H bc_fill.H)
    list(APPEND _sources DefineVelocity.cpp  face_velocity.H  Kernels.H  main.cpp  Tagging.H)
    list(APPEND _sources Src_K/Adv_K.H  Src_K/compute_flux_${D}D_K.H  Src_K/slope_K.H)
    list(TRANSFORM _sources PREPEND Source/)
    list(APPEND _sources Exec/Prob.H)

    # List of input files
    set(_input_files inputs-ci)
    list(TRANSFORM _input_files PREPEND "Exec/")

    setup_test(${D} _sources _input_files)

    unset( _sources )
    unset( _input_files   )
endforeach()
