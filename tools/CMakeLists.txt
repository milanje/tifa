#------------------------------------------------------------------------------
# Builds the TIFA command-line factorization programs and the parameter-tuning
# tool.
#------------------------------------------------------------------------------

set(TIFA_TOOLS_INCLUDE_DIRS
    "${PROJECT_BINARY_DIR}"                        # generated tifa_config.h (<tifa_config.h>)
    "${PROJECT_SOURCE_DIR}/lib/algo/include"
    "${PROJECT_SOURCE_DIR}/lib/data/include"
    "${PROJECT_SOURCE_DIR}/lib/utils/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

#------------------------------------------------------------------------------
# Sources shared by every program -> compiled once.
#------------------------------------------------------------------------------
add_library(tools_common OBJECT
    src/tool_utils.c
    src/factoring_program.c
    src/common_funcs.c
)
target_include_directories(tools_common PRIVATE ${TIFA_TOOLS_INCLUDE_DIRS})
target_link_libraries(tools_common PRIVATE GMP::GMP)   # for gmp.h at compile time

#------------------------------------------------------------------------------
# One executable per algorithm, and tifa_factor, which chooses them. Each is
# named tifa-<name> (tifa-ecm, tifa-p-minus-1, tifa-factor): plain names such
# as ecm clash with other tools.
#------------------------------------------------------------------------------
set(TIFA_PROGRAMS
    cfrac ecm fermat siqs squfof rho tdiv p_minus_1 p_plus_1 tifa_factor)
foreach(_prog IN LISTS TIFA_PROGRAMS)
    add_executable(${_prog}
        src/${_prog}_program.c
        $<TARGET_OBJECTS:tools_common>
    )
    target_include_directories(${_prog} PRIVATE ${TIFA_TOOLS_INCLUDE_DIRS})
    target_link_libraries(${_prog} PRIVATE Tifa::tifa_static)
    string(REGEX REPLACE "^tifa_" "" _name "${_prog}")
    string(REPLACE "_" "-" _name "${_name}")
    set_target_properties(${_prog} PROPERTIES OUTPUT_NAME "tifa-${_name}")
endforeach()

install(TARGETS ${TIFA_PROGRAMS}
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")

#------------------------------------------------------------------------------
# The parameter-tuning tool (tifa-tune).
#------------------------------------------------------------------------------
add_subdirectory(tune)
