if (NOT DEFINED TEST_EXECUTABLE OR NOT DEFINED TEST_INPUT
    OR NOT DEFINED EXPECTED_REGEX)
    message(FATAL_ERROR "Expected-failure test is missing a required argument")
endif ()

execute_process(
    COMMAND "${TEST_EXECUTABLE}" "${TEST_INPUT}"
    RESULT_VARIABLE test_result
    OUTPUT_VARIABLE test_stdout
    ERROR_VARIABLE test_stderr
)

if (test_result STREQUAL "0")
    message(FATAL_ERROR "Command unexpectedly succeeded: ${TEST_EXECUTABLE} ${TEST_INPUT}")
endif ()

set(test_output "${test_stdout}\n${test_stderr}")
if (NOT test_output MATCHES "${EXPECTED_REGEX}")
    message(FATAL_ERROR
        "Command failed for the wrong reason. Expected '${EXPECTED_REGEX}'.\n"
        "Result: ${test_result}\nOutput:\n${test_output}")
endif ()
