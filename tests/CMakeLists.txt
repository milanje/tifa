#------------------------------------------------------------------------------
# Test suite. Dependency-free: each test is a small C program whose exit code
# is the pass/fail signal, registered with CTest.
#------------------------------------------------------------------------------

set(TIFA_TEST_INCLUDE_DIRS
    "${PROJECT_BINARY_DIR}"                     # generated tifa_config.h
    "${PROJECT_SOURCE_DIR}/lib/algo/include"
    "${PROJECT_SOURCE_DIR}/lib/data/include"
    "${PROJECT_SOURCE_DIR}/lib/utils/include"
    "${CMAKE_CURRENT_SOURCE_DIR}"               # tifa_test.h
)

function(tifa_add_test name)
    add_executable(${name} ${name}.c)
    target_include_directories(${name} PRIVATE ${TIFA_TEST_INCLUDE_DIRS})
    target_link_libraries(${name} PRIVATE tifa_static)
    add_test(NAME ${name} COMMAND ${name})
    # Cap runtime so a regression that hangs (e.g. a bad loop) is reported as a
    # failure rather than stalling the run indefinitely.
    set_tests_properties(${name} PROPERTIES TIMEOUT 120)
endfunction()

tifa_add_test(test_funcs)
tifa_add_test(test_matrix)
tifa_add_test(test_tifa_factor)
tifa_add_test(test_algorithms)
tifa_add_test(test_siqs)
tifa_add_test(test_factor_u64)
tifa_add_test(test_mont_ecm)
tifa_add_test(test_mont_pm1)
tifa_add_test(test_cofactor)

# The tuning-loader test reads fixture files from the source tree; hand it their
# location as a compile definition (no env var, so the test is hermetic).
tifa_add_test(test_tuning)
target_compile_definitions(test_tuning PRIVATE
    TIFA_TEST_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures")

tifa_add_test(test_threads)
target_link_libraries(test_threads PRIVATE Threads::Threads)
target_compile_definitions(test_threads PRIVATE
    TIFA_TEST_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures")

# The tifa_factor command-line program: a 60-bit semiprime, and
# 2^3 . 3^2 . 1000000007^2 . 1000000009 (repeated factors, several sizes).
add_test(NAME tifa_factor_program_semiprime
         COMMAND tifa_factor 1000000016000000063)
set_tests_properties(tifa_factor_program_semiprime PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "1000000007[^0-9].*1000000009[^0-9].*Factorization is complete")
add_test(NAME tifa_factor_program_powers
         COMMAND tifa_factor 72000001656000012600000031752)
set_tests_properties(tifa_factor_program_powers PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "2 \\^ 3[^0-9].*3 \\^ 2[^0-9].*1000000009 \\^ 1[^0-9].*1000000007 \\^ 2[^0-9].*Factorization is complete")
