if (NOT DEFINED TEST_EXECUTABLE OR NOT DEFINED TEST_INPUT
    OR NOT DEFINED EXPECTED_REGEX)
    message(FATAL_ERROR "Expected-failure test is missing a required argument")
endif ()

set(test_command ${TEST_PREFIX} "${TEST_EXECUTABLE}" "${TEST_INPUT}" ${TEST_SUFFIX})

execute_process(
    COMMAND ${test_command}
    RESULT_VARIABLE test_result
    OUTPUT_VARIABLE test_stdout
    ERROR_VARIABLE test_stderr
)

if (test_result STREQUAL "0")
    message(FATAL_ERROR "Command unexpectedly succeeded: ${test_command}")
endif ()

set(test_output "${test_stdout}\n${test_stderr}")
if (NOT test_output MATCHES "${EXPECTED_REGEX}")
    message(FATAL_ERROR
        "Command failed for the wrong reason. Expected '${EXPECTED_REGEX}'.\n"
        "Result: ${test_result}\nOutput:\n${test_output}")
endif ()
