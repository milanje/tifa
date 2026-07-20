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
# GMP_FOUND, GMP_INCLUDE_DIRS and GMP_LIBRARIES.

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

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
    REQUIRED_VARS GMP_LIBRARY GMP_INCLUDE_DIR
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
