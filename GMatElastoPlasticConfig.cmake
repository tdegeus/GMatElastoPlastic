# GMatElastoPlastic cmake module
#
# This module sets the target:
#
#     GMatElastoPlastic
#
# In addition, it sets the following variables:
#
#     GMatElastoPlastic_FOUND - true if the library is found
#     GMatElastoPlastic_VERSION - the library's version
#     GMatElastoPlastic_INCLUDE_DIRS - directory containing the library's headers
#
# The following support targets are defined to simplify things:
#
#     GMatElastoPlastic::compiler_warnings - enable compiler warnings
#     GMatElastoPlastic::assert - enable library assertions
#     GMatElastoPlastic::debug - enable all assertions (slow)

include(CMakeFindDependencyMacro)

# Define target "GMatElastoPlastic"

if(NOT TARGET GMatElastoPlastic)
    include("${CMAKE_CURRENT_LIST_DIR}/GMatElastoPlasticTargets.cmake")
endif()

# Define "GMatElastoPlastic_INCLUDE_DIRS"

get_target_property(
    GMatElastoPlastic_INCLUDE_DIRS
    GMatElastoPlastic
    INTERFACE_INCLUDE_DIRECTORIES)

# Find dependencies

find_dependency(xtensor)

# Define support target "GMatElastoPlastic::compiler_warnings"

if(NOT TARGET GMatElastoPlastic::compiler_warnings)
    add_library(GMatElastoPlastic::compiler_warnings INTERFACE IMPORTED)
    if(MSVC)
        set_property(
            TARGET GMatElastoPlastic::compiler_warnings
            PROPERTY INTERFACE_COMPILE_OPTIONS
            /W4)
    else()
        set_property(
            TARGET GMatElastoPlastic::compiler_warnings
            PROPERTY INTERFACE_COMPILE_OPTIONS
            -Wall -Wextra -pedantic -Wno-unknown-pragmas)
    endif()
endif()

# Define support target "GMatElastoPlastic::assert"

if(NOT TARGET GMatElastoPlastic::assert)
    add_library(GMatElastoPlastic::assert INTERFACE IMPORTED)
    set_property(
        TARGET GMatElastoPlastic::assert
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        GMATELASTOPLASTIC_ENABLE_ASSERT)
endif()

# Define support target "GMatElastoPlastic::debug"

if(NOT TARGET GMatElastoPlastic::debug)
    add_library(GMatElastoPlastic::debug INTERFACE IMPORTED)
    set_property(
        TARGET GMatElastoPlastic::debug
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        XTENSOR_ENABLE_ASSERT GMATELASTOPLASTIC_ENABLE_ASSERT)
endif()
