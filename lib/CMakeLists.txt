#------------------------------------------------------------------------------
# Builds the TIFA library: a static (libtifa.a) and a shared (libtifa.so /
# libtifa.dylib) version from a single set of compiled objects.
#------------------------------------------------------------------------------

set(TIFA_SOURCES
    # Algorithms
    algo/src/cfrac.c
    algo/src/cofactor.c
    algo/src/ecm.c
    algo/src/fermat.c
    algo/src/siqs.c
    algo/src/squfof.c
    algo/src/rho.c
    algo/src/p_minus_1.c
    algo/src/p_plus_1.c
    algo/src/tdiv.c
    algo/src/tifa_factor.c
    algo/src/tuning.c
    # Precomputed data (generated tables)
    data/src/first_primes.c
    data/src/z2z8_curves.c
    # Structures, utilities, etc.
    utils/src/approx.c
    utils/src/array.c
    utils/src/bernsteinisms.c
    utils/src/factoring_machine.c
    utils/src/factor_u64.c
    utils/src/funcs.c
    utils/src/gauss_elim.c
    utils/src/gmp_utils.c
    utils/src/hashtable.c
    utils/src/lindep.c
    utils/src/linked_list.c
    utils/src/lp_table.c
    utils/src/matrix.c
    utils/src/ecm_schedule.c
    utils/src/mont_ecm.c
    utils/src/res_tdiv.c
    utils/src/siqs_poly.c
    utils/src/siqs_sieve.c
    utils/src/smooth_filter.c
    utils/src/sqrt_cont_frac.c
    utils/src/stopwatch.c
    utils/src/x_array_list.c
    utils/src/x_tree.c
)

# The baked tuning tables are compiled in only when explicitly enabled; the
# committed file is otherwise an inert placeholder (see tuning.h / tuning.c).
if(TIFA_USE_BAKED_TUNING)
    list(APPEND TIFA_SOURCES algo/src/tuning_baked.c)
endif()

# Include directories needed to compile the library sources.
set(TIFA_BUILD_INCLUDE_DIRS
    "${PROJECT_BINARY_DIR}"                        # generated tifa_config.h
    "${CMAKE_CURRENT_SOURCE_DIR}/algo/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/data/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils/include"
)

#------------------------------------------------------------------------------
# Compile the sources once into an object library.
#------------------------------------------------------------------------------
add_library(tifa_objects OBJECT ${TIFA_SOURCES})
# Only the declarations marked TIFA_API (the public headers) are exported by
# the shared library.
set_target_properties(tifa_objects PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    C_VISIBILITY_PRESET hidden
)
target_include_directories(tifa_objects PRIVATE ${TIFA_BUILD_INCLUDE_DIRS})
# GMP::GMP brings in gmp.h's include directory for compilation.
target_link_libraries(tifa_objects PRIVATE GMP::GMP)
if(TIFA_USE_GMP_INTERNAL_FUNCS AND GMP_INTERNAL_INCDIR)
    target_include_directories(tifa_objects PRIVATE "${GMP_INTERNAL_INCDIR}")
endif()

#------------------------------------------------------------------------------
# Wrap the objects into a static and a shared library.
#------------------------------------------------------------------------------
add_library(tifa_static STATIC $<TARGET_OBJECTS:tifa_objects>)
add_library(tifa_shared SHARED $<TARGET_OBJECTS:tifa_objects>)

foreach(_t tifa_static tifa_shared)
    set_target_properties(${_t} PROPERTIES OUTPUT_NAME tifa)
    target_include_directories(${_t} PUBLIC
        "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/algo/include>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/data/include>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/utils/include>"
    )
    target_link_libraries(${_t} PUBLIC GMP::GMP Threads::Threads)
    if(TIFA_M_LIBRARY)
        # By name: the found path would end up in the exported targets.
        target_link_libraries(${_t} PUBLIC m)
    endif()
endforeach()

# Below 1.0, a minor version may break the ABI (0.3 removed a field from
# cfrac_params_t, for one), so it is part of the soname.
if(PROJECT_VERSION_MAJOR EQUAL 0)
    set(TIFA_SOVERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")
else()
    set(TIFA_SOVERSION "${PROJECT_VERSION_MAJOR}")
endif()
set_target_properties(tifa_shared PROPERTIES
    VERSION "${PROJECT_VERSION}"
    SOVERSION "${TIFA_SOVERSION}"
)

# Exported (consumable) names: Tifa::tifa (shared), Tifa::tifa_static (static).
set_target_properties(tifa_shared PROPERTIES EXPORT_NAME tifa)
set_target_properties(tifa_static PROPERTIES EXPORT_NAME tifa_static)

# In-tree alias used by the tools/ subdirectory.
add_library(Tifa::tifa_static ALIAS tifa_static)

#------------------------------------------------------------------------------
# Public headers (installed flat into <prefix>/include/tifa)
#------------------------------------------------------------------------------
set(TIFA_PUBLIC_HEADERS
    # Algorithms
    algo/include/cfrac.h
    algo/include/cofactor.h
    algo/include/ecm.h
    algo/include/fermat.h
    algo/include/siqs.h
    algo/include/squfof.h
    algo/include/rho.h
    algo/include/p_minus_1.h
    algo/include/p_plus_1.h
    algo/include/tdiv.h
    algo/include/tifa_factor.h
    algo/include/tuning.h
    # Data
    data/include/first_primes.h
    # Structures and utilities
    utils/include/array.h
    utils/include/bitstring_t.h
    utils/include/exit_codes.h
    utils/include/factoring_machine.h
    utils/include/factor_u64.h
    utils/include/methods.h
    # Umbrella header
    tifa.h
)

install(FILES ${TIFA_PUBLIC_HEADERS}
    DESTINATION "${TIFA_INSTALL_INCLUDEDIR}")
install(FILES "${PROJECT_BINARY_DIR}/tifa_config.h"
    DESTINATION "${TIFA_INSTALL_INCLUDEDIR}")
