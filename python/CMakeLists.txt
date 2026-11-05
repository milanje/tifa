#------------------------------------------------------------------------------
# The extension module tifa._tifa, built by pip through scikit-build-core
# (pyproject.toml at the top of the tree) against the limited C API of Python
# 3.10, so that one wheel serves every later version.
#------------------------------------------------------------------------------
if(CMAKE_VERSION VERSION_LESS 3.26)
    message(FATAL_ERROR "TIFA_BUILD_PYTHON needs CMake 3.26 or later")
endif()

find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.SABIModule)

Python_add_library(_tifa MODULE USE_SABI 3.10 WITH_SOABI _tifa.c)
target_link_libraries(_tifa PRIVATE tifa_static)
# A static GMP (the wheels') keeps its symbols to the module.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_link_options(_tifa PRIVATE "LINKER:--exclude-libs,ALL")
endif()

install(TARGETS _tifa LIBRARY DESTINATION tifa COMPONENT python)
if(DEFINED SKBUILD_SCRIPTS_DIR)
    install(TARGETS tifa_tune RUNTIME DESTINATION "${SKBUILD_SCRIPTS_DIR}"
            COMPONENT python)
endif()
