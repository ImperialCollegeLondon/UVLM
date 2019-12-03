cmake_minimum_required(VERSION 3.9)
project(uvlm)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP)

file(GLOB SOURCES src/cpp_interface.cpp)
include_directories(include)
add_library(uvlm SHARED ${SOURCES})

set_property(TARGET uvlm PROPERTY CXX_STANDARD 17)

if(OpenMP_CXX_FOUND)
    target_link_libraries(uvlm PUBLIC OpenMP::OpenMP_CXX Eigen3::Eigen)
else()
    target_link_libraries(uvlm Eigen3::Eigen)
endif()

# Add custom debug flags
target_compile_options(uvlm PUBLIC $<$<CONFIG:DEBUG>: -Wall>)

# For all compilers
if(CMAKE_BUILD_TYPE STREQUAL Release)
    find_package(BLAS)
endif()
target_compile_options(uvlm PUBLIC $<$<CONFIG:RELEASE>: -fomit-frame-pointer -ffast-math -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE>)

install(TARGETS uvlm
    LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)

