
cmake_minimum_required(VERSION 3.0)

# Basic settings
# ==============

project(GMatElastoPlastic)

option(BUILD_TESTS "${PROJECT_NAME} Build tests" OFF)

# Version
# =======

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/GMatElastoPlastic/config.h" _GMAT_ver
     REGEX "#define GMATELASTOPLASTIC_VERSION_(MAJOR|MINOR|PATCH)")

foreach(ver ${_GMAT_ver})
    if(ver MATCHES "#define GMATELASTOPLASTIC_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
        set(_GMAT_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
    endif()
endforeach()

set(_GMAT_VERSION ${_GMAT_VERSION_MAJOR}.${_GMAT_VERSION_MINOR}.${_GMAT_VERSION_PATCH})

message(STATUS "Building GMatElastoPlastic v${_GMAT_VERSION}")

# Set target
# ==========

find_package(xtensor REQUIRED)
find_package(GMatTensor REQUIRED)
find_package(GMatElastic REQUIRED)

add_library(GMatElastoPlastic INTERFACE)

target_include_directories(GMatElastoPlastic INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_link_libraries(GMatElastoPlastic INTERFACE xtensor GMatTensor GMatElastic)

# Installation
# ============

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/" DESTINATION include)

install(TARGETS GMatElastoPlastic EXPORT GMatElastoPlastic-targets)

install(
    EXPORT GMatElastoPlastic-targets
    FILE GMatElastoPlasticTargets.cmake
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GMatElastoPlastic")

set(_GMAT ${CMAKE_SIZEOF_VOID_P})
unset(CMAKE_SIZEOF_VOID_P)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/GMatElastoPlasticConfigVersion.cmake"
    VERSION ${_GMAT_VERSION}
    COMPATIBILITY AnyNewerVersion)

set(CMAKE_SIZEOF_VOID_P ${_GMAT})

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/GMatElastoPlasticConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/GMatElastoPlasticConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GMatElastoPlastic")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GMatElastoPlastic.pc.in"
               "${CMAKE_CURRENT_BINARY_DIR}/GMatElastoPlastic.pc" @ONLY)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/GMatElastoPlastic.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")

# Add builds
# ==========

include("GMatElastoPlasticConfig.cmake")

set(CTEST_TEST_TARGET_ALIAS RUN_ALL_TESTS)
include(CTest)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
