cmake_minimum_required(VERSION 2.6)
project(colour_cube C CXX)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release (-DCMAKE_BUILD_TYPE=...)")
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(NO_REFLECTION ON)
set(NO_FFTW ON)
set(NO_GSL ON)
set(NO_LAPACK ON)
set(NO_WCSLIB ON)

find_package(vif REQUIRED)

# Include dependencies
message(STATUS ${VIF_INCLUDE_DIRS})
include_directories(${VIF_INCLUDE_DIRS})
include_directories(${COMMON_INCLUDE_DIR})

# Build ColourCube
add_executable(colour_cube colour_cube.cpp)
target_link_libraries(colour_cube ${VIF_LIBRARIES})
install(TARGETS colour_cube DESTINATION bin)

# Build ColourCube cube builder
add_executable(make_colour_cube make_colour_cube.cpp)
target_link_libraries(make_colour_cube ${VIF_LIBRARIES})
install(TARGETS make_colour_cube DESTINATION bin)
