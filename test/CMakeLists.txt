cmake_minimum_required(VERSION 3.0)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    project(GMatElastoPlastic-test)
    find_package(GMatElastoPlastic REQUIRED CONFIG)
endif()

option(WARNINGS "Show build warnings" ON)
option(ASSERT "Enable assertions" ON)
option(XSIMD "Use xsimd and 'march=native' optimisations" OFF)

set(CMAKE_BUILD_TYPE Release)

find_package(Catch2 REQUIRED)
find_package(xtensor REQUIRED)
add_executable(main main.cpp Cartesian3d.cpp)
target_link_libraries(main PRIVATE Catch2::Catch2 GMatElastoPlastic)

if(WARNINGS)
    target_link_libraries(main PRIVATE GMatElastoPlastic::compiler_warnings)
endif()

if(ASSERT)
    target_link_libraries(main PRIVATE GMatElastoPlastic::assert)
endif()

if(XSIMD)
    target_link_libraries(main PRIVATE xtensor::optimize xtensor::use_xsimd)
endif()
