if (NOT (3 IN_LIST AMReX_SPACEDIM))
    return()
endif ()

set(_sources  main.cpp)

set(_input_files
    inputs.cube.mc
    inputs.cube.legacy
    inputs.cubes.mc
    inputs.cubes.legacy
    inputs.bunny
    inputs.bunny.legacy
    inputs.multilevel.mc
    inputs.multilevel.rebuild.mc
    inputs.nodal.mc
    inputs.open-stl.mc
    cube.stl
    cubes.stl
    stanford_bunny.stl
    STANFORD_BUNNY_LICENSE.md
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

add_mc_test_variant(EB_MarchingCubesMultilevel_3d inputs.multilevel.mc)
add_mc_test_variant(EB_MarchingCubesMultilevelRebuild_3d inputs.multilevel.rebuild.mc)
add_mc_test_variant(EB_MarchingCubesNodal_3d inputs.nodal.mc)
add_mc_test_variant(EB_MarchingCubesCubeLegacy_3d inputs.cube.legacy)
add_mc_test_variant(EB_MarchingCubesCubes_3d inputs.cubes.mc)
add_mc_test_variant(EB_MarchingCubesCubesLegacy_3d inputs.cubes.legacy)
add_mc_test_variant(EB_MarchingCubesBunny_3d inputs.bunny)
add_mc_test_variant(EB_MarchingCubesBunnyLegacy_3d inputs.bunny.legacy)

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
    inputs.open-stl.mc "STL is open, nonmanifold, or inconsistently oriented")

unset(_sources)
unset(_input_files)
unset(_mc_test_labels)
