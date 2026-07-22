#------------------------------------------------------------------------------
# TIFA - Tools for Integer FActorization
#
# CMake build system (replaces the historical SCons build).
#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.20)

project(Tifa
    VERSION 0.1.5
    DESCRIPTION "TIFA - Tools for Integer FActorization"
    LANGUAGES C
)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CheckTypeSize)
include(CMakePackageConfigHelpers)

#------------------------------------------------------------------------------
# Build type / language standard
#------------------------------------------------------------------------------
# Default to an optimized build, matching the historical SCons "-O3" default.
# Note: we deliberately keep the Release flags free of -DNDEBUG for parity with
# the old build (which never defined NDEBUG).
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
        STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()
set(CMAKE_C_FLAGS_RELEASE "-O3")

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

if(NOT DEFINED CMAKE_POSITION_INDEPENDENT_CODE)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Install headers into a dedicated tifa/ subfolder (the sources include their
# headers by bare basename, so consumers add <prefix>/include/tifa to -I).
set(TIFA_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/tifa")

#------------------------------------------------------------------------------
# Options (successors of BuildOptions.py)
#------------------------------------------------------------------------------
option(TIFA_USE_GMP_INTERNAL_FUNCS "Use some of GMP's internal functions (needs gmp-impl.h)" OFF)
option(TIFA_USE_OWN_RAND           "Use TIFA's basic pseudo-random generator"                OFF)
option(TIFA_USE_CALLOC_MEMSET      "Use calloc/memset for zero-init of integer arrays"       ON)
option(TIFA_PRINT_ERROR            "Print critical error messages on stderr"                 ON)
option(TIFA_ALLOW_VERBOSE          "Master switch enabling per-algorithm verbose output"     OFF)
option(TIFA_ALLOW_TIMING           "Master switch enabling per-algorithm timing messages"    OFF)
option(TIFA_ALLOW_EXTRA_TIMING     "Enable extra (potentially costly) timing measurements"   OFF)

set(TIFA_ALGOS CFRAC ECM FERMAT SIQS SQUFOF RHO TDIV PM1 PP1)
foreach(_a IN LISTS TIFA_ALGOS)
    string(TOLOWER "${_a}" _al)
    option(TIFA_VERBOSE_${_a} "Verbose output for the ${_al} function (needs TIFA_ALLOW_VERBOSE)" ON)
    option(TIFA_TIMING_${_a}  "Timing messages for the ${_al} function (needs TIFA_ALLOW_TIMING)"  ON)
endforeach()

set(TIFA_BITSTRING_T "uint64_t" CACHE STRING "Native C type used to represent a string of bits")

set(GMP_INCDIR          "" CACHE PATH "Directory containing gmp.h (optional hint for FindGMP)")
set(GMP_LIBDIR          "" CACHE PATH "Directory containing the GMP library (optional hint)")
set(GMP_INTERNAL_INCDIR "" CACHE PATH "Directory containing gmp-impl.h (used only with TIFA_USE_GMP_INTERNAL_FUNCS)")

# Documentation options
option(TIFA_GENERATE_DOXYGEN_DOC  "Provide a 'doc' target that runs Doxygen" ON)
option(TIFA_GENERATE_DOXYGEN_HTML "Doxygen: generate HTML output"            ON)
option(TIFA_GENERATE_DOXYGEN_PDF  "Doxygen: generate PDF (LaTeX) output"     ON)
option(TIFA_GENERATE_DOXYGEN_XML  "Doxygen: generate XML output"             OFF)
option(TIFA_GENERATE_DOXYGEN_RTF  "Doxygen: generate RTF output"             OFF)
option(TIFA_GENERATE_MAN_PAGES    "Generate Perl module man pages via pod2man" ON)
set(TIFA_DOXYGEN_PAPER_TYPE "a4wide" CACHE STRING "Doxygen PDF paper type")
set_property(CACHE TIFA_DOXYGEN_PAPER_TYPE PROPERTY
    STRINGS a4 a4wide letter legal executive)

option(TIFA_BUILD_TESTS "Build the test suite (CTest)" ON)
option(TIFA_WERROR "Treat compiler warnings as errors (used in CI)" OFF)

#------------------------------------------------------------------------------
# Dependencies
#------------------------------------------------------------------------------
find_package(GMP REQUIRED)
find_library(TIFA_M_LIBRARY m)   # libm; may be absent (folded into libc) on some platforms

#------------------------------------------------------------------------------
# Feature detection (compile-only; cross-compilation safe, no run-probes)
#------------------------------------------------------------------------------
set(CMAKE_EXTRA_INCLUDE_FILES stdint.h stdlib.h)
check_type_size("${TIFA_BITSTRING_T}"  TIFA_CFG_SIZEOF_BITSTRING_T LANGUAGE C)
check_type_size("unsigned long int"    TIFA_CFG_SIZEOF_ULONG_T     LANGUAGE C)
check_type_size("long int"             TIFA_SIZEOF_LONG            LANGUAGE C)
unset(CMAKE_EXTRA_INCLUDE_FILES)

if(NOT TIFA_CFG_SIZEOF_BITSTRING_T)
    message(FATAL_ERROR
        "Could not determine sizeof(${TIFA_BITSTRING_T}). "
        "Set -DTIFA_BITSTRING_T=<a valid C integer type>.")
endif()

# Precomputed integer square roots of LONG_MAX / ULONG_MAX (used in runtime
# comparisons). Only 32- and 64-bit longs are supported; extend if needed.
if(TIFA_SIZEOF_LONG EQUAL 8)
    set(TIFA_CFG_SQRT_LONG_MAX "3037000499L")
elseif(TIFA_SIZEOF_LONG EQUAL 4)
    set(TIFA_CFG_SQRT_LONG_MAX "46340L")
else()
    message(FATAL_ERROR "Unsupported sizeof(long int) = ${TIFA_SIZEOF_LONG}")
endif()
if(TIFA_CFG_SIZEOF_ULONG_T EQUAL 8)
    set(TIFA_CFG_SQRT_ULONG_MAX "4294967295UL")
elseif(TIFA_CFG_SIZEOF_ULONG_T EQUAL 4)
    set(TIFA_CFG_SQRT_ULONG_MAX "65535UL")
else()
    message(FATAL_ERROR "Unsupported sizeof(unsigned long int) = ${TIFA_CFG_SIZEOF_ULONG_T}")
endif()

#------------------------------------------------------------------------------
# Assemble the substitution values for tifa_config.h
#------------------------------------------------------------------------------
# Helper: map the option named by flag_var to the literal 1 or 0.
function(tifa_bool01 out_var flag_var)
    if(${flag_var})
        set(${out_var} 1 PARENT_SCOPE)
    else()
        set(${out_var} 0 PARENT_SCOPE)
    endif()
endfunction()

set(TIFA_CFG_SHORTNAME "TIFA")
set(TIFA_CFG_FULLNAME  "TIFA - Tools for Integer FActorization")
set(TIFA_CFG_TARNAME   "tifa")
set(TIFA_CFG_STRING    "TIFA v${PROJECT_VERSION}")
set(TIFA_CFG_BUGREPORT "milanje at gmail dot com")

set(TIFA_CFG_VERSION       "${PROJECT_VERSION}")
set(TIFA_CFG_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
set(TIFA_CFG_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(TIFA_CFG_VERSION_PATCH "${PROJECT_VERSION_PATCH}")

math(EXPR _tifa_wordsize_bits "${TIFA_CFG_SIZEOF_ULONG_T} * 8")
set(TIFA_CFG_COMPILED_ON_OS      "${CMAKE_HOST_SYSTEM_NAME}")
set(TIFA_CFG_COMPILED_ON_RELEASE "${CMAKE_HOST_SYSTEM_VERSION}")
set(TIFA_CFG_COMPILED_ON_MACHINE "${CMAKE_HOST_SYSTEM_PROCESSOR}")
set(TIFA_CFG_COMPILED_ON
    "${CMAKE_HOST_SYSTEM_NAME}-${CMAKE_HOST_SYSTEM_VERSION}-${CMAKE_HOST_SYSTEM_PROCESSOR}-${_tifa_wordsize_bits}bit")

tifa_bool01(TIFA_CFG_USE_CALLOC_MEMSET      TIFA_USE_CALLOC_MEMSET)
tifa_bool01(TIFA_CFG_USE_GMP_INTERNAL_FUNCS TIFA_USE_GMP_INTERNAL_FUNCS)
tifa_bool01(TIFA_CFG_USE_OWN_RAND           TIFA_USE_OWN_RAND)
tifa_bool01(TIFA_CFG_PRINT_ERROR            TIFA_PRINT_ERROR)

set(TIFA_CFG_BITSTRING_T "${TIFA_BITSTRING_T}")

# Per-algorithm verbose/timing, gated by the master switches (matches SConstruct).
foreach(_a IN LISTS TIFA_ALGOS)
    if(TIFA_ALLOW_VERBOSE AND TIFA_VERBOSE_${_a})
        set(TIFA_CFG_VERBOSE_${_a} 1)
    else()
        set(TIFA_CFG_VERBOSE_${_a} 0)
    endif()
    if(TIFA_ALLOW_TIMING AND TIFA_TIMING_${_a})
        set(TIFA_CFG_TIMING_${_a} 1)
    else()
        set(TIFA_CFG_TIMING_${_a} 0)
    endif()
endforeach()
if(TIFA_ALLOW_TIMING AND TIFA_ALLOW_EXTRA_TIMING)
    set(TIFA_CFG_ALLOW_EXTRA_TIMING 1)
else()
    set(TIFA_CFG_ALLOW_EXTRA_TIMING 0)
endif()

configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/tifa_config.h.in"
    "${PROJECT_BINARY_DIR}/tifa_config.h"
    @ONLY
)

#------------------------------------------------------------------------------
# Warning flags for TIFA's own compilation. Applied here (after all configure-
# time feature checks) so they affect only the targets built below, never the
# check_type_size()/find_package() probes, and are not exported to consumers.
#------------------------------------------------------------------------------
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra)
    if(TIFA_WERROR)
        add_compile_options(-Werror)
    endif()
endif()

#------------------------------------------------------------------------------
# Build the library and the tools
#------------------------------------------------------------------------------
add_subdirectory(lib)
add_subdirectory(tools)

#------------------------------------------------------------------------------
# Tests (CTest)
#------------------------------------------------------------------------------
if(TIFA_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

#------------------------------------------------------------------------------
# Installation & packaging
#------------------------------------------------------------------------------
install(TARGETS tifa_static tifa_shared
    EXPORT TifaTargets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    INCLUDES DESTINATION "${TIFA_INSTALL_INCLUDEDIR}"
)

install(EXPORT TifaTargets
    FILE TifaTargets.cmake
    NAMESPACE Tifa::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Tifa"
)

configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/TifaConfig.cmake.in"
    "${PROJECT_BINARY_DIR}/TifaConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Tifa"
)
write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/TifaConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion
)
install(FILES
    "${PROJECT_BINARY_DIR}/TifaConfig.cmake"
    "${PROJECT_BINARY_DIR}/TifaConfigVersion.cmake"
    "${PROJECT_SOURCE_DIR}/cmake/FindGMP.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Tifa"
)

# pkg-config. TIFA's public headers include <gmp.h>, so the GMP dependency
# must be expressed here. Prefer delegating to GMP's own gmp.pc when one was
# found at configure time; otherwise embed the discovered GMP paths so the
# .pc is self-contained.
if(PC_GMP_FOUND)
    set(TIFA_PC_REQUIRES_PRIVATE "gmp")
    set(TIFA_PC_LIBS_PRIVATE "-lm")
    set(TIFA_PC_CFLAGS_EXTRA "")
else()
    set(TIFA_PC_REQUIRES_PRIVATE "")
    get_filename_component(_tifa_gmp_libdir "${GMP_LIBRARY}" DIRECTORY)
    set(TIFA_PC_LIBS_PRIVATE "-L${_tifa_gmp_libdir} -lgmp -lm")
    set(TIFA_PC_CFLAGS_EXTRA " -I${GMP_INCLUDE_DIR}")
endif()
configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/tifa.pc.in"
    "${PROJECT_BINARY_DIR}/tifa.pc"
    @ONLY
)
install(FILES "${PROJECT_BINARY_DIR}/tifa.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

#------------------------------------------------------------------------------
# Documentation (Doxygen) - replaces the old os.system() based doc builder
#------------------------------------------------------------------------------
if(TIFA_GENERATE_DOXYGEN_DOC)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        set(_doc_dir "${PROJECT_BINARY_DIR}/doc")
        file(MAKE_DIRECTORY "${_doc_dir}")

        macro(_tifa_yesno out flag)
            if(${flag})
                set(${out} "YES")
            else()
                set(${out} "NO")
            endif()
        endmacro()
        _tifa_yesno(_doxy_html  TIFA_GENERATE_DOXYGEN_HTML)
        _tifa_yesno(_doxy_latex TIFA_GENERATE_DOXYGEN_PDF)
        _tifa_yesno(_doxy_xml   TIFA_GENERATE_DOXYGEN_XML)
        _tifa_yesno(_doxy_rtf   TIFA_GENERATE_DOXYGEN_RTF)

        # The historical Doxyfile template uses %TOKEN% placeholders; substitute
        # them at configure time (all values are known now).
        file(READ "${PROJECT_SOURCE_DIR}/doxygen.conf.in" _doxy)
        string(REPLACE "%PROJECT%"        "${TIFA_CFG_SHORTNAME}"      _doxy "${_doxy}")
        string(REPLACE "%VERSION%"        "${PROJECT_VERSION}"         _doxy "${_doxy}")
        string(REPLACE "%DOCDIR%"         "${_doc_dir}"                _doxy "${_doxy}")
        string(REPLACE "%TOPSRCDIR%"      "${PROJECT_SOURCE_DIR}"      _doxy "${_doxy}")
        string(REPLACE "%PAPER_TYPE%"     "${TIFA_DOXYGEN_PAPER_TYPE}" _doxy "${_doxy}")
        string(REPLACE "%GENERATE_HTML%"  "${_doxy_html}"              _doxy "${_doxy}")
        string(REPLACE "%GENERATE_LATEX%" "${_doxy_latex}"             _doxy "${_doxy}")
        string(REPLACE "%GENERATE_XML%"   "${_doxy_xml}"               _doxy "${_doxy}")
        string(REPLACE "%GENERATE_RTF%"   "${_doxy_rtf}"               _doxy "${_doxy}")
        file(WRITE "${_doc_dir}/doxygen.conf" "${_doxy}")

        add_custom_target(doc
            COMMAND "${DOXYGEN_EXECUTABLE}" "${_doc_dir}/doxygen.conf"
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM
        )

        if(TIFA_GENERATE_DOXYGEN_PDF)
            find_package(LATEX COMPONENTS PDFLATEX QUIET)
            if(LATEX_FOUND)
                # Doxygen emits a LaTeX Makefile; build the PDF best-effort.
                add_custom_command(TARGET doc POST_BUILD
                    COMMAND sh -c "if [ -d '${_doc_dir}/latex' ]; then make -C '${_doc_dir}/latex' || echo 'TIFA: PDF build failed (non-fatal).'; fi"
                    VERBATIM
                )
            else()
                message(STATUS "TIFA: LaTeX not found; the 'doc' target will not build a PDF.")
            endif()
        endif()
    else()
        message(STATUS "TIFA: Doxygen not found; the 'doc' target is unavailable.")
    endif()
endif()

#------------------------------------------------------------------------------
# Source packaging (replaces the broken SCons 'dist' target)
#   cpack --config CPackSourceConfig.cmake   (or: cmake --build . --target package_source)
#------------------------------------------------------------------------------
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "tifa-${PROJECT_VERSION}")
set(CPACK_VERBATIM_VARIABLES YES)
set(CPACK_SOURCE_IGNORE_FILES
    "/\\.git/"
    "/build/"
    "/\\.sconf_temp/"
    "\\.sconsign"
    "\\.sconsign\\.dblite"
    "\\.o$"
    "\\.a$"
    "\\.so$"
    "\\.dylib$"
    "\\.tgz$"
    "\\.tar\\.gz$"
    "/tifa_config\\.h$"
)
include(CPack)

#------------------------------------------------------------------------------
# Configuration summary
#------------------------------------------------------------------------------
message(STATUS "")
message(STATUS "TIFA ${PROJECT_VERSION} configured:")
message(STATUS "  Build type          : ${CMAKE_BUILD_TYPE}")
message(STATUS "  BITSTRING_T         : ${TIFA_BITSTRING_T} (sizeof=${TIFA_CFG_SIZEOF_BITSTRING_T})")
message(STATUS "  GMP include dir     : ${GMP_INCLUDE_DIRS}")
message(STATUS "  GMP library         : ${GMP_LIBRARIES}")
message(STATUS "  libm                : ${TIFA_M_LIBRARY}")
message(STATUS "  Install prefix      : ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
