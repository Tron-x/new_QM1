cmake_minimum_required(VERSION 2.8)
project(qm)

set(CXX_FLAGS
    -g
    -std=c++11
    -O3
    )

string(REPLACE ";" " " CMAKE_CXX_FLAGS "${CXX_FLAGS}")

add_subdirectory(src)
