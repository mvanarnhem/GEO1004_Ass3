cmake_minimum_required(VERSION 3.26)
project(src)

set(CMAKE_CXX_STANDARD 17)
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

add_executable(src main.cpp)
