if (NOT (3 IN_LIST AMReX_SPACEDIM))
    return()
endif ()

set(_sources  main.cpp)

set(_input_files
    inputs
    inputs.cleanup
    inputs.bunny
    inputs.bunny.legacy
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
    stanford_bunny.stl
    STANFORD_BUNNY_LICENSE.md
)

setup_test(3 _sources _input_files)

add_test_variant(EB_MarchingCubesCleanup_3d
    Test_EB_MarchingCubes_3d inputs.cleanup)
add_test_variant(EB_MarchingCubesBunnyMC_3d
    Test_EB_MarchingCubes_3d inputs.bunny)
add_test_variant(EB_MarchingCubesBunnyLegacy_3d
    Test_EB_MarchingCubes_3d inputs.bunny.legacy)
add_test_variant(EB_MarchingCubesMultilevel_3d
    Test_EB_MarchingCubes_3d inputs.multilevel)
add_test_variant(EB_MarchingCubesMultilevelRebuild_3d
    Test_EB_MarchingCubes_3d inputs.multilevel.rebuild)
add_test_variant(EB_MarchingCubesNodal_3d
    Test_EB_MarchingCubes_3d inputs.nodal)
add_test_variant(EB_MarchingCubesPeriodic_3d
    Test_EB_MarchingCubes_3d inputs.periodic)
add_test_variant(EB_MarchingCubesScaled_3d
    Test_EB_MarchingCubes_3d inputs.scaled)
add_test_variant(EB_MarchingCubesExtendDomain_3d
    Test_EB_MarchingCubes_3d inputs.extend-domain)
add_test_variant(EB_MarchingCubesExtendDomainLegacy_3d
    Test_EB_MarchingCubes_3d inputs.extend-domain
    eb2.stl_geometry_method=legacy)

function(add_expected_failure _test_name _inputs_filename _expected_regex)
    add_test(
        NAME ${_test_name}
        COMMAND ${CMAKE_COMMAND}
            -DTEST_EXECUTABLE=$<TARGET_FILE:Test_EB_MarchingCubes_3d>
            -DTEST_INPUT=${_inputs_filename}
            -DEXPECTED_REGEX=${_expected_regex}
            -P ${CMAKE_CURRENT_LIST_DIR}/expect_failure.cmake
        WORKING_DIRECTORY $<TARGET_FILE_DIR:Test_EB_MarchingCubes_3d>
    )
endfunction()

add_expected_failure(EB_MarchingCubesRejectAnisotropic_3d
    inputs.anisotropic "requires a 3D Cartesian grid")
add_expected_failure(EB_MarchingCubesRejectOpenSTL_3d
    inputs.open-stl "STL is open, nonmanifold, or inconsistently oriented")

unset(_sources)
unset(_input_files)
