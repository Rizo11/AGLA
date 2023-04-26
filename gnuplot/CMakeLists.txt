cmake_minimum_required(VERSION 3.25)
project(gnuplot)

set(CMAKE_CXX_STANDARD 17)

add_executable(gnuplot main.cpp)

find_package(Gnuplot REQUIRED)
include_directories(${GNUPLOT_INCLUDE_DIRS})
target_link_libraries(gnuplot ${GNUPLOT_LIBRARIES})