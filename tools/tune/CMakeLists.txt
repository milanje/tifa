#------------------------------------------------------------------------------
# tifa_tune: the parameter-tuning tool. Links the static library and calls its
# algorithms directly; each run is forked under a watchdog unless --timeout 0.
#------------------------------------------------------------------------------
add_executable(tifa_tune
    tune.c
    tune_grid.c
    tune_numgen.c
    tune_algos.c
)
target_include_directories(tifa_tune PRIVATE
    ${TIFA_TOOLS_INCLUDE_DIRS}
    "${CMAKE_CURRENT_SOURCE_DIR}"
)
target_link_libraries(tifa_tune PRIVATE Tifa::tifa_static GMP::GMP)
if(TIFA_M_LIBRARY)
    target_link_libraries(tifa_tune PRIVATE ${TIFA_M_LIBRARY})   # libm for fmax/floor
endif()

install(TARGETS tifa_tune
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")

#------------------------------------------------------------------------------
# Convenience "quick'n'dirty" tuning pass into <build>/tuning. Not part of ALL.
#------------------------------------------------------------------------------
add_custom_target(tune-quick
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${CMAKE_BINARY_DIR}/tuning"
    COMMAND tifa_tune --quick --out "${CMAKE_BINARY_DIR}/tuning"
    DEPENDS tifa_tune
    COMMENT "Running a quick tuning pass into ${CMAKE_BINARY_DIR}/tuning"
    VERBATIM)

#------------------------------------------------------------------------------
# Regenerate the committed baked tuning tables (a full tuning pass, then bake).
# Opt-in because it can be lengthy; the developer commits the regenerated files.
#------------------------------------------------------------------------------
if(TIFA_REGEN_TUNING)
    add_custom_target(regen-tuning
        COMMAND "${CMAKE_COMMAND}" -E make_directory "${CMAKE_BINARY_DIR}/tuning"
        COMMAND tifa_tune --out "${CMAKE_BINARY_DIR}/tuning"
        COMMAND tifa_tune
                --emit-c "${PROJECT_SOURCE_DIR}/lib/algo/src/tuning_baked.c"
                --emit-h "${PROJECT_SOURCE_DIR}/lib/algo/include/tuning_baked.h"
                --in "${CMAKE_BINARY_DIR}/tuning"
        DEPENDS tifa_tune
        COMMENT "Regenerating baked tuning tables (full tuning pass)"
        VERBATIM)
endif()
