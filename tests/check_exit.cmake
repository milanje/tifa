#------------------------------------------------------------------------------
# Runs PROGRAM with the ;-separated ARGS (and an empty standard input) and
# fails unless it exits with EXPECTED.
#------------------------------------------------------------------------------
execute_process(
    COMMAND "${PROGRAM}" ${ARGS}
    INPUT_FILE /dev/null
    OUTPUT_QUIET
    ERROR_QUIET
    RESULT_VARIABLE status
)
if(NOT status STREQUAL EXPECTED)
    message(FATAL_ERROR "${PROGRAM} ${ARGS}: exit status ${status}, expected ${EXPECTED}")
endif()
