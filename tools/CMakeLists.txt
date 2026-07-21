#------------------------------------------------------------------------------
# Builds the TIFA command-line factorization programs and installs the Perl
# tooling (benchmark/plot scripts, modules, man pages).
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
# One executable per algorithm.
#------------------------------------------------------------------------------
set(TIFA_PROGRAMS cfrac ecm fermat siqs squfof rho tdiv p_minus_1)
foreach(_prog IN LISTS TIFA_PROGRAMS)
    add_executable(${_prog}
        src/${_prog}_program.c
        $<TARGET_OBJECTS:tools_common>
    )
    target_include_directories(${_prog} PRIVATE ${TIFA_TOOLS_INCLUDE_DIRS})
    target_link_libraries(${_prog} PRIVATE Tifa::tifa_static)
endforeach()

install(TARGETS ${TIFA_PROGRAMS}
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")

#------------------------------------------------------------------------------
# Perl scripts, configuration files and modules
#------------------------------------------------------------------------------
install(PROGRAMS
    scripts/factorize.pl
    scripts/benchmarker.pl
    scripts/plotmaker.pl
    scripts/mergeres.pl
    scripts/extractres.pl
    scripts/gencomposites.pl
    DESTINATION "${CMAKE_INSTALL_BINDIR}")

install(FILES
    scripts/factorize.conf
    scripts/benchmarker.conf
    scripts/plotmaker.conf
    DESTINATION "${CMAKE_INSTALL_DATADIR}/tifa/conf")

set(TIFA_INSTALL_PERLMODDIR "${CMAKE_INSTALL_DATADIR}/tifa/perl"
    CACHE PATH "Where to install TIFA's Perl modules")

set(TIFA_PERL_MODULES
    Bencher
    DataDescriptor
    FormatConverter
    GnuPlotter
    NumberGenerator
    Program
    ProgramRepository
    Rotor
    SimpleConfigReader
)

set(_pm_files)
foreach(_m IN LISTS TIFA_PERL_MODULES)
    list(APPEND _pm_files "scripts/Tifa/${_m}.pm")
endforeach()
install(FILES ${_pm_files}
    DESTINATION "${TIFA_INSTALL_PERLMODDIR}/Tifa")

#------------------------------------------------------------------------------
# Perl module man pages (via pod2man, if available)
#------------------------------------------------------------------------------
find_program(POD2MAN_EXECUTABLE pod2man)
if(TIFA_GENERATE_MAN_PAGES AND POD2MAN_EXECUTABLE)
    set(_man_files)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/man")
    foreach(_m IN LISTS TIFA_PERL_MODULES)
        set(_src "${CMAKE_CURRENT_SOURCE_DIR}/scripts/Tifa/${_m}.pm")
        # Avoid ':' in build-tree paths (a colon breaks the Makefile generator);
        # the installed page keeps the conventional Tifa::<Module>.3pm name.
        set(_out "${CMAKE_CURRENT_BINARY_DIR}/man/Tifa-${_m}.3pm")
        add_custom_command(
            OUTPUT "${_out}"
            COMMAND "${POD2MAN_EXECUTABLE}"
                    -c "TIFA Documentation"
                    -r "TIFA v${PROJECT_VERSION}"
                    -n "Tifa::${_m}"
                    "${_src}" "${_out}"
            DEPENDS "${_src}"
            COMMENT "Generating man page Tifa::${_m}.3pm"
            VERBATIM
        )
        list(APPEND _man_files "${_out}")
        install(FILES "${_out}"
            RENAME "Tifa::${_m}.3pm"
            DESTINATION "${CMAKE_INSTALL_MANDIR}/man3")
    endforeach()
    add_custom_target(man ALL DEPENDS ${_man_files})
elseif(TIFA_GENERATE_MAN_PAGES)
    message(STATUS "TIFA: pod2man not found; skipping Perl module man pages.")
endif()
