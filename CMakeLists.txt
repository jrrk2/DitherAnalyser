cmake_minimum_required(VERSION 3.16)
project(DitherAnalyser VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt6 REQUIRED COMPONENTS Widgets Charts Concurrent)
find_package(PkgConfig REQUIRED)
pkg_check_modules(CFITSIO REQUIRED cfitsio)

# Use the system Qt6 build of stellarsolver (not the Qt5 one in /usr/local)
set(STELLARSOLVER_INCLUDE_DIRS /usr/include/libstellarsolver)
set(STELLARSOLVER_LIBRARIES /usr/lib/x86_64-linux-gnu/libstellarsolver6.so)

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
    ${STELLARSOLVER_LIBRARIES}
)
