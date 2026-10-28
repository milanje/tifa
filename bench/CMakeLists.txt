#------------------------------------------------------------------------------
# The benchmark programs, built with -DTIFA_BUILD_BENCH=ON and never installed:
# siqs-bench and factor-bench time the library, siqs-check checks SIQS on
# numbers of random sizes, and ab.sh compares two builds (see README.md).
#------------------------------------------------------------------------------
foreach(_prog siqs_bench factor_bench siqs_check)
    add_executable(${_prog} ${_prog}.c)
    target_include_directories(${_prog} PRIVATE
        "${PROJECT_BINARY_DIR}"
        "${PROJECT_SOURCE_DIR}/lib/algo/include"
        "${PROJECT_SOURCE_DIR}/lib/data/include"
        "${PROJECT_SOURCE_DIR}/lib/utils/include"
    )
    target_link_libraries(${_prog} PRIVATE Tifa::tifa_static GMP::GMP)
    string(REPLACE "_" "-" _name "${_prog}")
    set_target_properties(${_prog} PROPERTIES OUTPUT_NAME "${_name}")
endforeach()

if(TIFA_BUILD_TESTS)
    add_test(NAME siqs_check COMMAND siqs_check 50 120 60)
endif()
