function(run_checked)
    execute_process(
        COMMAND ${ARGN}
        RESULT_VARIABLE _result
        OUTPUT_VARIABLE _stdout
        ERROR_VARIABLE _stderr
    )

    if(NOT _result EQUAL 0)
        string(JOIN " " _command ${ARGN})
        message(FATAL_ERROR
            "Command failed (${_result}): ${_command}\n"
            "stdout:\n${_stdout}\n"
            "stderr:\n${_stderr}")
    endif()
endfunction()

function(assert_exists _path)
    if(NOT EXISTS "${_path}" AND NOT IS_SYMLINK "${_path}")
        message(FATAL_ERROR "Expected installed file is missing: ${_path}")
    endif()
endfunction()

file(REMOVE_RECURSE
    "${AMREX_TEST_INSTALL_PREFIX}"
    "${AMREX_TEST_INSTALL_STAGE_ROOT}"
    "${AMREX_TEST_INSTALL_TEST_BUILD_DIR}")
file(MAKE_DIRECTORY
    "${AMREX_TEST_INSTALL_PREFIX}"
    "${AMREX_TEST_INSTALL_STAGE_ROOT}")

set(_install_cmd
    "${CMAKE_COMMAND}" "--install" "${AMREX_TEST_INSTALL_BUILD_DIR}"
    "--prefix" "${AMREX_TEST_INSTALL_PREFIX}")
if(AMREX_TEST_INSTALL_CONFIG)
    list(APPEND _install_cmd "--config" "${AMREX_TEST_INSTALL_CONFIG}")
endif()
run_checked(${_install_cmd})

set(_legacy_lib "${AMREX_TEST_INSTALL_PREFIX}/lib/${AMREX_TEST_INSTALL_LEGACY_LIB_NAME}")
assert_exists("${_legacy_lib}")

if(CMAKE_HOST_WIN32)
    message(STATUS "Skipping DESTDIR staged install check on Windows.")
else()
    set(_stage_install_cmd
        "${CMAKE_COMMAND}" "-E" "env" "DESTDIR=${AMREX_TEST_INSTALL_STAGE_ROOT}"
        "${CMAKE_COMMAND}" "--install" "${AMREX_TEST_INSTALL_BUILD_DIR}"
        "--prefix" "${AMREX_TEST_INSTALL_STAGE_PREFIX}")
    if(AMREX_TEST_INSTALL_CONFIG)
        list(APPEND _stage_install_cmd "--config" "${AMREX_TEST_INSTALL_CONFIG}")
    endif()
    run_checked(${_stage_install_cmd})

    set(_staged_legacy_lib
        "${AMREX_TEST_INSTALL_STAGE_ROOT}${AMREX_TEST_INSTALL_STAGE_PREFIX}/lib/${AMREX_TEST_INSTALL_LEGACY_LIB_NAME}")
    assert_exists("${_staged_legacy_lib}")

    file(READ "${AMREX_TEST_INSTALL_BUILD_DIR}/install_manifest.txt" _manifest)
    set(_manifest_legacy_lib
        "${AMREX_TEST_INSTALL_STAGE_PREFIX}/lib/${AMREX_TEST_INSTALL_LEGACY_LIB_NAME}")
    string(FIND "${_manifest}" "${_manifest_legacy_lib}" _manifest_legacy_lib_pos)
    if(_manifest_legacy_lib_pos EQUAL -1)
        message(FATAL_ERROR
            "Install manifest is missing the legacy library entry: ${_manifest_legacy_lib}")
    endif()

    set(_staged_manifest_legacy_lib
        "${AMREX_TEST_INSTALL_STAGE_ROOT}${_manifest_legacy_lib}")
    string(FIND "${_manifest}" "${_staged_manifest_legacy_lib}" _staged_manifest_legacy_lib_pos)
    if(NOT _staged_manifest_legacy_lib_pos EQUAL -1)
        message(FATAL_ERROR
            "Install manifest should not record DESTDIR paths: ${_staged_manifest_legacy_lib}")
    endif()
endif()

set(_configure_cmd
    "${CMAKE_COMMAND}"
    "-S" "${AMREX_TEST_INSTALL_TEST_PROJECT_DIR}"
    "-B" "${AMREX_TEST_INSTALL_TEST_BUILD_DIR}"
    "-DAMReX_ROOT=${AMREX_TEST_INSTALL_PREFIX}"
    "-DCMAKE_CXX_COMPILER=${AMREX_TEST_INSTALL_CXX_COMPILER}")
if(AMREX_TEST_INSTALL_FORTRAN_COMPILER)
    list(APPEND _configure_cmd
        "-DCMAKE_Fortran_COMPILER=${AMREX_TEST_INSTALL_FORTRAN_COMPILER}")
endif()
if(NOT AMREX_TEST_INSTALL_IS_MULTI_CONFIG AND AMREX_TEST_INSTALL_CONFIG)
    list(APPEND _configure_cmd
        "-DCMAKE_BUILD_TYPE=${AMREX_TEST_INSTALL_CONFIG}")
endif()
run_checked(${_configure_cmd})

set(_build_cmd
    "${CMAKE_COMMAND}" "--build" "${AMREX_TEST_INSTALL_TEST_BUILD_DIR}")
if(AMREX_TEST_INSTALL_CONFIG)
    list(APPEND _build_cmd "--config" "${AMREX_TEST_INSTALL_CONFIG}")
endif()
run_checked(${_build_cmd})

file(REMOVE_RECURSE "${AMREX_TEST_INSTALL_TEST_BUILD_DIR}")
