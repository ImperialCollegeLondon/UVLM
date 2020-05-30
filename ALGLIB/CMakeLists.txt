cmake_minimum_required(VERSION 3.9)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/Modules/")
project(alglib)

#if(NOT CMAKE_BUILD_TYPE)
#  set(CMAKE_BUILD_TYPE Release)
#endif()

file(GLOB SOURCES src/*.cpp)
add_library(alglib ${SOURCES})
target_include_directories (alglib PUBLIC src)

set_property(TARGET alglib PROPERTY CXX_STANDARD 14)

if(CMAKE_BUILD_TYPE STREQUAL Release)
    find_package(BLAS)
    set(BUILD_SHARED_LIBS ON)
    find_package(MKL)
    if(NOT MKL_FOUND)
        find_package(LAPACK REQUIRED)
        target_link_libraries(alglib PUBLIC ${LAPACK_LIBRARIES})
    else()
        include_directories(${MKL_INCLUDE_DIR})
        target_link_libraries(alglib PUBLIC ${MKL_LIBRARIES})
    endif()
endif()

find_package(Eigen3 REQUIRED)
target_link_libraries(alglib PUBLIC Eigen3::Eigen)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(alglib PUBLIC OpenMP::OpenMP_CXX)
endif()

find_package(Eigen3 REQUIRED)
target_link_libraries(alglib PUBLIC Eigen3::Eigen)

# For all compilers
target_compile_options(alglib PUBLIC $<$<CONFIG:RELEASE>: -fomit-frame-pointer -ffast-math> -fPIC)

# Dodgy thing for MacOS Catalina with anyting other than CLANG++
# Apple has removed the /usr/include dir in the folder structure.
# So compilers that are not the original Clang (AppleClang in CMAKE)
# do not work without this line
if(APPLE)
    if(NOT ${CMAKE_SYSTEM_VERSION} VERSION_LESS "19.0.0")
       if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
            message(FATAL_ERROR "Intel icpc compiler with MacOS Catalina does not work properly as I am writing this.  If you think this is wrong (or outdated), please file an issue in github.com/imperialcollegelondon/sharpy.  For now, try another C++ compiler, like g++ or clang")
       endif()
        if(NOT ${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
            target_compile_options(alglib PUBLIC -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk)
        endif()
    endif()
endif()

# Add custom debug flags
#target_compile_options(alglib PUBLIC $<$<CONFIG:DEBUG>: -Wall>)

install(TARGETS alglib
   LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib
   ARCHIVE DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)

