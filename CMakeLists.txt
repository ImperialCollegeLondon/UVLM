cmake_minimum_required(VERSION 3.9)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")
project(uvlm)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

file(GLOB SOURCES src/cpp_interface.cpp)
include_directories(include)
add_library(uvlm SHARED ${SOURCES})

set_property(TARGET uvlm PROPERTY CXX_STANDARD 17)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(uvlm PUBLIC OpenMP::OpenMP_CXX)
endif()

find_package(Eigen3 REQUIRED)
target_link_libraries(uvlm PUBLIC Eigen3::Eigen)

# For all compilers
if(CMAKE_BUILD_TYPE STREQUAL Release)
    SET(BUILD_SHARED_LIBS ON)
    find_package(MKL REQUIRED)
    include_directories(${MKL_INCLUDE_DIR})
    target_link_libraries(uvlm PUBLIC ${MKL_LIBRARIES})
endif()
target_compile_options(uvlm PUBLIC $<$<CONFIG:RELEASE>: -fomit-frame-pointer -ffast-math -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE>)

# Add custom debug flags
target_compile_options(uvlm PUBLIC $<$<CONFIG:DEBUG>: -Wall>)

install(TARGETS uvlm
    LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)

