cmake_minimum_required(VERSION 2.8.9)
project(emd-suite)

include_directories(include)

file(GLOB SOURCES "src/*.cpp" "src/*.c")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_executable(emd-suite ${SOURCES})
