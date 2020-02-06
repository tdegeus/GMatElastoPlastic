cmake_minimum_required(VERSION 3.0)

# Basic settings
# ==============

project(GMatElastoPlastic)

option(BUILD_TESTS "${PROJECT_NAME} Build tests" OFF)

# Version
# =======

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/GMatElastoPlastic/config.h" GMatElastoPlastic_version_defines
     REGEX "#define GMATELASTOPLASTIC_VERSION_(MAJOR|MINOR|PATCH)")

foreach(ver ${GMatElastoPlastic_version_defines})
    if(ver MATCHES "#define GMATELASTOPLASTIC_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
    set(GMATELASTOPLASTIC_VERSION_${CMAKE_MATCH_1}
        "${CMAKE_MATCH_2}"
        CACHE INTERNAL "")
    endif()
endforeach()

set(GMATELASTOPLASTIC_VERSION
    ${GMATELASTOPLASTIC_VERSION_MAJOR}.${GMATELASTOPLASTIC_VERSION_MINOR}.${GMATELASTOPLASTIC_VERSION_PATCH})

message(STATUS "Building GMatElastoPlastic v${GMATELASTOPLASTIC_VERSION}")

# Set target
# ==========

find_package(xtensor REQUIRED)

add_library(GMatElastoPlastic INTERFACE)

target_include_directories(GMatElastoPlastic INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_link_libraries(GMatElastoPlastic INTERFACE xtensor)

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

set(_GMATELASTOPLASTIC ${CMAKE_SIZEOF_VOID_P})
unset(CMAKE_SIZEOF_VOID_P)

write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/GMatElastoPlasticConfigVersion.cmake"
    VERSION ${GMATELASTOPLASTIC_VERSION}
    COMPATIBILITY AnyNewerVersion)

set(CMAKE_SIZEOF_VOID_P ${_GMATELASTOPLASTIC})

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

if(BUILD_TESTS)
    add_subdirectory(test)
endif()
