cmake_minimum_required(VERSION 3.16)
project(DitherAnalyser VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt6 REQUIRED COMPONENTS Widgets Charts)
find_package(PkgConfig REQUIRED)
pkg_check_modules(CFITSIO REQUIRED cfitsio)

add_executable(DitherAnalyser
    main.cpp
    DitherAnalyser.cpp
    DitherAnalyser.h
)

target_include_directories(DitherAnalyser PRIVATE ${CFITSIO_INCLUDE_DIRS})
target_link_libraries(DitherAnalyser PRIVATE
    Qt6::Widgets
    Qt6::Charts
    ${CFITSIO_LIBRARIES}
)
