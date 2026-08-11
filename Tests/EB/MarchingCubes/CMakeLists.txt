if (NOT (3 IN_LIST AMReX_SPACEDIM))
    return()
endif ()

set(_sources  main.cpp)

set(_input_files
    inputs
    inputs.cleanup
    inputs.multilevel
    inputs.multilevel.rebuild
    inputs.nodal
    inputs.periodic
    inputs.scaled
    inputs.extend-domain
    inputs.anisotropic
    inputs.open-stl
    cube.stl
    open_surface.stl
)

setup_test(3 _sources _input_files COMBINE_MPI_OMP NTASKS 2 NTHREADS 2)

set(_mc_test_labels EB MarchingCubes)
if (NOT AMReX_GPU_BACKEND STREQUAL NONE)
    list(APPEND _mc_test_labels GPU)
endif()
set_tests_properties(EB_MarchingCubes_3d PROPERTIES
    LABELS "${_mc_test_labels}")

function(add_mc_test_variant _test_name _inputs_filename)
    set(_cmd $<TARGET_FILE:Test_EB_MarchingCubes_3d> ${_inputs_filename}
        algorithm_tests=0 ${ARGN})
    if (AMReX_MPI)
        set(_cmd ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            ${MPIEXEC_PREFLAGS} ${_cmd} ${MPIEXEC_POSTFLAGS})
    endif()
    add_test(
        NAME ${_test_name}
        COMMAND ${_cmd}
        WORKING_DIRECTORY $<TARGET_FILE_DIR:Test_EB_MarchingCubes_3d>
    )
    if (AMReX_OMP)
        set_tests_properties(${_test_name} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=2)
    endif()
    set_tests_properties(${_test_name} PROPERTIES LABELS "${_mc_test_labels}")
endfunction()

add_mc_test_variant(EB_MarchingCubesCleanup_3d inputs.cleanup)
add_mc_test_variant(EB_MarchingCubesMultilevel_3d inputs.multilevel)
add_mc_test_variant(EB_MarchingCubesMultilevelRebuild_3d inputs.multilevel.rebuild)
add_mc_test_variant(EB_MarchingCubesNodal_3d inputs.nodal)
add_mc_test_variant(EB_MarchingCubesPeriodic_3d inputs.periodic)
add_mc_test_variant(EB_MarchingCubesScaled_3d inputs.scaled)
add_mc_test_variant(EB_MarchingCubesExtendDomain_3d inputs.extend-domain)
add_mc_test_variant(EB_MarchingCubesExtendDomainLegacy_3d inputs.extend-domain
    eb2.stl_geometry_method=legacy)
add_mc_test_variant(EB_MarchingCubesLegacySTL_3d inputs
    narrow_band_test=0
    eb2.stl_geometry_method=legacy
    eb2.eb_surface_stl_file=cube_legacy.stl)

find_package(Python3 COMPONENTS Interpreter REQUIRED)
add_test(
    NAME EB_MarchingCubesCompareSTL_3d
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/compare_stl.py
        cube_legacy.stl cube_mc.stl
        --relative-tolerance 0.0001
        --absolute-tolerance 0.00001
    WORKING_DIRECTORY $<TARGET_FILE_DIR:Test_EB_MarchingCubes_3d>
)
set_tests_properties(
    EB_MarchingCubes_3d
    EB_MarchingCubesLegacySTL_3d
    PROPERTIES FIXTURES_SETUP MarchingCubesSTL
)
set_tests_properties(EB_MarchingCubesCompareSTL_3d PROPERTIES
    FIXTURES_REQUIRED MarchingCubesSTL
    LABELS "${_mc_test_labels}"
)

function(add_expected_failure _test_name _inputs_filename _expected_regex)
    set(_test_prefix)
    set(_test_suffix)
    if (AMReX_MPI)
        set(_test_prefix ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            ${MPIEXEC_PREFLAGS})
        set(_test_suffix ${MPIEXEC_POSTFLAGS})
    endif()
    add_test(
        NAME ${_test_name}
        COMMAND ${CMAKE_COMMAND}
            "-DTEST_PREFIX=${_test_prefix}"
            "-DTEST_SUFFIX=${_test_suffix}"
            -DTEST_EXECUTABLE=$<TARGET_FILE:Test_EB_MarchingCubes_3d>
            -DTEST_INPUT=${_inputs_filename}
            -DEXPECTED_REGEX=${_expected_regex}
            -P ${CMAKE_CURRENT_LIST_DIR}/expect_failure.cmake
        WORKING_DIRECTORY $<TARGET_FILE_DIR:Test_EB_MarchingCubes_3d>
    )
    if (AMReX_OMP)
        set_tests_properties(${_test_name} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=2)
    endif()
    set_tests_properties(${_test_name} PROPERTIES LABELS "${_mc_test_labels}")
endfunction()

add_expected_failure(EB_MarchingCubesRejectOpenSTL_3d
    inputs.open-stl "STL is open, nonmanifold, or inconsistently oriented")
add_expected_failure(EB_MarchingCubesRejectAnisotropic_3d
    inputs.anisotropic "requires a 3D Cartesian grid")

unset(_sources)
unset(_input_files)
unset(_mc_test_labels)
