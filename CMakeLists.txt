cmake_minimum_required(VERSION 3.9)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")
project(uvlm)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Eigen3 REQUIRED)
find_package(OpenMP)

file(GLOB SOURCES src/cpp_interface.cpp)
include_directories(include)
add_library(uvlm SHARED ${SOURCES})

set_property(TARGET uvlm PROPERTY CXX_STANDARD 14)
target_link_libraries(uvlm PUBLIC Eigen3::Eigen)

if(OpenMP_CXX_FOUND)
    target_link_libraries(uvlm PUBLIC OpenMP::OpenMP_CXX)
endif()

if(CMAKE_BUILD_TYPE STREQUAL Release)
    find_package(BLAS)
    set(BUILD_SHARED_LIBS ON)
    find_package(MKL)
    if(NOT MKL_FOUND)
        find_package(LAPACK REQUIRED)
        target_link_libraries(uvlm PUBLIC ${LAPACK_LIBRARIES})
    else()
        include_directories(${MKL_INCLUDE_DIR})
        target_link_libraries(uvlm PUBLIC ${MKL_LIBRARIES})
    endif()
endif()

# Add custom debug flags
target_compile_options(uvlm PUBLIC $<$<CONFIG:DEBUG>: -Wall>)

# For all compilers
target_compile_options(uvlm PUBLIC $<$<CONFIG:RELEASE>: -fomit-frame-pointer -ffast-math -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE>)

install(TARGETS uvlm
    LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)

