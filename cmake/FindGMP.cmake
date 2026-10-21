# FindGMP.cmake - Locate the GNU Multiple Precision Arithmetic library (GMP).
#
# GMP does not reliably ship a pkg-config file, so this module searches the
# usual locations directly. It honours two optional hint variables kept from
# TIFA's historical build options:
#
#   GMP_INCDIR - directory containing gmp.h
#   GMP_LIBDIR - directory containing the GMP library
#
# On success it defines the imported target GMP::GMP as well as the variables
# GMP_FOUND, GMP_INCLUDE_DIRS, GMP_LIBRARIES and GMP_VERSION (read from gmp.h,
# and checked against the version asked for, if any).

if(TARGET GMP::GMP)
    set(GMP_FOUND TRUE)
    return()
endif()

# Use pkg-config only as a source of extra search hints (never required).
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(PC_GMP QUIET gmp)
endif()

find_path(GMP_INCLUDE_DIR
    NAMES gmp.h
    HINTS
        ${GMP_INCDIR}
        ${PC_GMP_INCLUDEDIR}
        ${PC_GMP_INCLUDE_DIRS}
    PATHS
        /opt/homebrew/include
        /usr/local/include
        /usr/include
)

find_library(GMP_LIBRARY
    NAMES gmp libgmp
    HINTS
        ${GMP_LIBDIR}
        ${PC_GMP_LIBDIR}
        ${PC_GMP_LIBRARY_DIRS}
    PATHS
        /opt/homebrew/lib
        /usr/local/lib
        /usr/lib
)

if(GMP_INCLUDE_DIR AND EXISTS "${GMP_INCLUDE_DIR}/gmp.h")
    file(STRINGS "${GMP_INCLUDE_DIR}/gmp.h" _gmp_version_lines
         REGEX "^#define[ \t]+__GNU_MP_VERSION(_MINOR|_PATCHLEVEL)?[ \t]+[0-9]+")
    foreach(_part "" "_MINOR" "_PATCHLEVEL")
        string(REGEX MATCH "__GNU_MP_VERSION${_part}[ \t]+([0-9]+)" _m
               "${_gmp_version_lines}")
        set(_gmp_v${_part} "${CMAKE_MATCH_1}")
    endforeach()
    if(NOT _gmp_v STREQUAL "")
        set(GMP_VERSION "${_gmp_v}.${_gmp_v_MINOR}.${_gmp_v_PATCHLEVEL}")
    endif()
    unset(_gmp_version_lines)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
    REQUIRED_VARS GMP_LIBRARY GMP_INCLUDE_DIR
    VERSION_VAR   GMP_VERSION
)

if(GMP_FOUND)
    set(GMP_LIBRARIES ${GMP_LIBRARY})
    set(GMP_INCLUDE_DIRS ${GMP_INCLUDE_DIR})
    if(NOT TARGET GMP::GMP)
        add_library(GMP::GMP UNKNOWN IMPORTED)
        set_target_properties(GMP::GMP PROPERTIES
            IMPORTED_LOCATION "${GMP_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}"
        )
    endif()
endif()

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY)
