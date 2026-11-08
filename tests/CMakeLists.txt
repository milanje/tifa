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
tifa_add_test(test_containers)
tifa_add_test(test_tifa_factor)
tifa_add_test(test_algorithms)
tifa_add_test(test_siqs)
tifa_add_test(test_factor_u64)
tifa_add_test(test_mont_ecm)
tifa_add_test(test_mont_pm1)
tifa_add_test(test_cofactor)
tifa_add_test(test_inputs)
tifa_add_test(test_prime_cert)
tifa_add_test(test_hostile_input)
target_link_libraries(test_hostile_input PRIVATE Threads::Threads)
tifa_add_test(test_time_limit)
target_link_libraries(test_time_limit PRIVATE Threads::Threads)
tifa_add_test(test_batch_gcd)
tifa_add_test(test_block_lanczos)

# The tuning-loader test reads fixture files from the source tree; hand it their
# location as a compile definition (no env var needed: the test unsets
# $TIFA_TUNING_DIR and never falls back to the installed tuning directory).
tifa_add_test(test_tuning)
target_compile_definitions(test_tuning PRIVATE
    TIFA_TEST_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures")

tifa_add_test(test_threads)
target_link_libraries(test_threads PRIVATE Threads::Threads)
target_compile_definitions(test_threads PRIVATE
    TIFA_TEST_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures")

# The tifa-factor command-line program: a 60-bit semiprime, and
# 2^3 . 3^2 . 1000000007^2 . 1000000009 (repeated factors, several sizes). The
# large factors' order depends on the build (with or without __int128).
add_test(NAME tifa_factor_program_semiprime
         COMMAND tifa_factor 1000000016000000063)
set_tests_properties(tifa_factor_program_semiprime PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "(1000000007[^0-9].*1000000009|1000000009[^0-9].*1000000007)[^0-9].*Factorization is complete")
add_test(NAME tifa_factor_program_powers
         COMMAND tifa_factor 72000001656000012600000031752)
set_tests_properties(tifa_factor_program_powers PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "2 \\^ 3[^0-9].*3 \\^ 2[^0-9].*(1000000009 \\^ 1[^0-9].*1000000007 \\^ 2|1000000007 \\^ 2[^0-9].*1000000009 \\^ 1)[^0-9].*Factorization is complete")
# A square of a composite: the program factors its root, then squares.
add_test(NAME tifa_factor_program_square
         COMMAND tifa_factor 1000000032000000382000002016000003969)
set_tests_properties(tifa_factor_program_square PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "(1000000007 \\^ 2[^0-9].*1000000009 \\^ 2|1000000009 \\^ 2[^0-9].*1000000007 \\^ 2)[^0-9].*Factorization is complete")
# In 4 threads: 576460752303423619 . 1152921504606847009 (119 bits), by SIQS.
add_test(NAME tifa_factor_program_threads
         COMMAND tifa_factor --threads 4 balanced 664613997892458106507825459650105571)
set_tests_properties(tifa_factor_program_threads PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "(576460752303423619[^0-9].*1152921504606847009|1152921504606847009[^0-9].*576460752303423619)[^0-9].*Factorization is complete")
# A hint: 340730379077795843 . 604614875560711445197 (128 bits).
add_test(NAME tifa_factor_program_hint
         COMMAND tifa_factor balanced 206010655745875572166197184221848916071)
set_tests_properties(tifa_factor_program_hint PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "hint: balanced.*(340730379077795843[^0-9].*604614875560711445197|604614875560711445197[^0-9].*340730379077795843)[^0-9].*Factorization is complete")
# --prove: a certificate with a step for the factor above 2^64, none for the
# other, then for a prime, 2^127 - 1.
add_test(NAME tifa_factor_program_prove
         COMMAND tifa_factor --prove balanced 206010655745875572166197184221848916071)
set_tests_properties(tifa_factor_program_prove PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "Primality certificates:.*(340730379077795843: \\[\\]\n.*604614875560711445197:\n\t    \\[\\[604614875560711445197, \\[\\[2, |604614875560711445197:\n\t    \\[\\[604614875560711445197, \\[\\[2, .*340730379077795843: \\[\\]\n).*Factorization is complete, its factors proven prime")
add_test(NAME tifa_factor_program_prove_prime
         COMMAND tifa_factor 170141183460469231731687303715884105727 --prove)
set_tests_properties(tifa_factor_program_prove_prime PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "170141183460469231731687303715884105727:\n\t    \\[\\[170141183460469231731687303715884105727, \\[\\[2, .*\n170141183460469231731687303715884105727 is prime\\.")
# The tifa-batch-gcd program: positions and gcds of the integers sharing a
# prime, every gcd with --all, hexadecimal with --hex.
set(_fixtures "${CMAKE_CURRENT_SOURCE_DIR}/fixtures")
add_test(NAME batch_gcd_program
         COMMAND batch_gcd "${_fixtures}/batch_gcd.txt")
set_tests_properties(batch_gcd_program PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "1 3\n2 3\n3 11\n4 11\n6 899\n7 899\n")
add_test(NAME batch_gcd_program_all
         COMMAND batch_gcd --all "${_fixtures}/batch_gcd.txt")
set_tests_properties(batch_gcd_program_all PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "3\n3\n11\n11\n1\n899\n899\n")
add_test(NAME batch_gcd_program_hex
         COMMAND batch_gcd --hex "${_fixtures}/batch_gcd_hex.txt")
set_tests_properties(batch_gcd_program_hex PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "1 3\n2 3\n3 b\n4 b\n")
# The algorithms' programs give the multiplicities of what they split.
add_test(NAME rho_program_square
         COMMAND rho 1000000032000000382000002016000003969)
set_tests_properties(rho_program_square PROPERTIES
    TIMEOUT 60
    PASS_REGULAR_EXPRESSION "1000000007 \\^ 2[^0-9].*1000000009 \\^ 2[^0-9].*Factorization is complete")
# A pass regex alone ignores the exit status: fail on sanitizer reports too.
set_tests_properties(tifa_factor_program_semiprime tifa_factor_program_powers
                     tifa_factor_program_square tifa_factor_program_hint
                     tifa_factor_program_prove tifa_factor_program_prove_prime
                     rho_program_square batch_gcd_program batch_gcd_program_all
                     batch_gcd_program_hex
                     PROPERTIES
    FAIL_REGULAR_EXPRESSION "ERROR|Sanitizer|runtime error")

# The programs' exit statuses (tool_utils.h): 0 complete (or a prime, or 1),
# 1 partial, 2 no factor found, 3 invalid input.
function(tifa_add_exit_test name expected program)
    add_test(NAME ${name}
             COMMAND ${CMAKE_COMMAND} "-DPROGRAM=$<TARGET_FILE:${program}>"
                     "-DARGS=${ARGN}" -DEXPECTED=${expected}
                     -P "${CMAKE_CURRENT_SOURCE_DIR}/check_exit.cmake")
    set_tests_properties(${name} PROPERTIES TIMEOUT 60)
endfunction()
tifa_add_exit_test(exit_complete  0 tifa_factor 1000000016000000063)
tifa_add_exit_test(exit_prime     0 tifa_factor 1000003)
tifa_add_exit_test(exit_partial   1 tdiv 10 3000000048000000189)
tifa_add_exit_test(exit_rho_complete 0 rho 1000000016000000063)
tifa_add_exit_test(exit_no_factor 2 tdiv 10 1000000016000000063)
tifa_add_exit_test(exit_zero      3 tifa_factor 0)
tifa_add_exit_test(exit_not_digit 3 tifa_factor 12x)
tifa_add_exit_test(exit_bad_argc  3 tifa_factor random 12 34)
tifa_add_exit_test(exit_bad_hint  3 tifa_factor 12 34)
tifa_add_exit_test(exit_no_input  3 tifa_factor)
tifa_add_exit_test(exit_help      0 tifa_factor --help)
tifa_add_exit_test(exit_prove_complete 0 tifa_factor --prove 1000000016000000063)
tifa_add_exit_test(exit_threads_complete 0 tifa_factor --threads 4 1000000016000000063)
tifa_add_exit_test(exit_threads_zero     3 tifa_factor --threads 0 1000000016000000063)
tifa_add_exit_test(exit_threads_missing  3 tifa_factor 1000000016000000063 --threads)
tifa_add_exit_test(exit_prove_prime    0 tifa_factor --prove 1000003)
tifa_add_exit_test(exit_prove_bad_hint 3 tifa_factor 12 --prove 34)
tifa_add_exit_test(exit_prove_help     0 tifa_factor --prove --help)
# Long numbers: 3 . (2^1279 - 1), 386 digits; 5000 digits are too many.
tifa_add_exit_test(exit_long_number    0 tifa_factor 31223796583993197245775720982092256615845786741800114415957337051210824179020080894036718325196972554527264841862843023374150698744970394724419719817031962288772388528167661037319550048984750651418878261643069433635081031442645436600887504153314780875197046655933482485615249755084766813998674495672756405541204666609421563080095306466198874086634152609208560702144162131667109506187261)
string(REPEAT "7" 5000 _too_long)
tifa_add_exit_test(exit_too_long       3 tifa_factor ${_too_long})
tifa_add_exit_test(exit_batch_gcd_shared  0 batch_gcd "${_fixtures}/batch_gcd.txt")
tifa_add_exit_test(exit_batch_gcd_none    2 batch_gcd "${_fixtures}/batch_gcd_coprime.txt")
tifa_add_exit_test(exit_batch_gcd_empty   2 batch_gcd)
tifa_add_exit_test(exit_batch_gcd_bad     3 batch_gcd "${_fixtures}/batch_gcd_bad.txt")
tifa_add_exit_test(exit_batch_gcd_missing 3 batch_gcd "${_fixtures}/no_such_file")
tifa_add_exit_test(exit_batch_gcd_option  3 batch_gcd --decimal)
tifa_add_exit_test(exit_batch_gcd_help    0 batch_gcd --help)

add_test(NAME too_many_digits COMMAND tifa_factor ${_too_long})
set_tests_properties(too_many_digits PROPERTIES
    TIMEOUT 60 PASS_REGULAR_EXPRESSION "more than 4095 digits")

# --version: "<program> (TIFA) <version>".
foreach(_prog tifa_factor siqs batch_gcd tifa_tune)
    add_test(NAME version_${_prog} COMMAND ${_prog} --version)
    set_tests_properties(version_${_prog} PROPERTIES
        TIMEOUT 60
        PASS_REGULAR_EXPRESSION "^tifa-[a-z-]+ \\(TIFA\\) ${PROJECT_VERSION}\n$")
endforeach()
tifa_add_exit_test(exit_help_h    0 siqs -h)
