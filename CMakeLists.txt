cmake_minimum_required(VERSION 3.16)
project(DitherAnalyser VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt6 REQUIRED COMPONENTS Widgets Charts Concurrent Network)
find_package(PkgConfig REQUIRED)
pkg_check_modules(CFITSIO REQUIRED cfitsio)

# Find stellarsolver — try cmake config first, fall back to manual search
# Prefer system path over /usr/local (which may have a Qt5 build)
find_package(StellarSolver QUIET
    PATHS /usr/lib/x86_64-linux-gnu/cmake /usr/lib/cmake
    NO_DEFAULT_PATH)
if(NOT StellarSolver_FOUND AND NOT TARGET StellarSolver::stellarsolver)
    find_package(StellarSolver QUIET)
endif()
if(StellarSolver_FOUND OR TARGET StellarSolver::stellarsolver)
    set(STELLARSOLVER_TARGET StellarSolver::stellarsolver)
else()
    # Fallback for system installs without cmake config
    find_path(STELLARSOLVER_INCLUDE_DIR stellarsolver.h
        PATHS
            /usr/include/libstellarsolver
            /usr/local/include/libstellarsolver
        PATH_SUFFIXES libstellarsolver
    )
    find_library(STELLARSOLVER_LIBRARY
        NAMES stellarsolver6 stellarsolver
        PATHS /usr/lib/x86_64-linux-gnu /usr/local/lib
    )
    if(STELLARSOLVER_INCLUDE_DIR AND STELLARSOLVER_LIBRARY)
        set(STELLARSOLVER_INCLUDE_DIRS ${STELLARSOLVER_INCLUDE_DIR})
        set(STELLARSOLVER_TARGET ${STELLARSOLVER_LIBRARY})
    else()
        message(FATAL_ERROR "stellarsolver not found")
    endif()
endif()

add_executable(DitherAnalyser
    main.cpp
    DitherAnalyser.cpp
    DitherAnalyser.h
)

target_include_directories(DitherAnalyser PRIVATE
    ${CFITSIO_INCLUDE_DIRS}
    ${STELLARSOLVER_INCLUDE_DIRS}
)
target_link_libraries(DitherAnalyser PRIVATE
    Qt6::Widgets
    Qt6::Charts
    Qt6::Concurrent
    ${CFITSIO_LIBRARIES}
    ${STELLARSOLVER_TARGET}
)
