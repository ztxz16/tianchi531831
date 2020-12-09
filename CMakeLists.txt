cmake_minimum_required(VERSION 3.5)

project(test LANGUAGES CXX)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -O2")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(opt opt.cpp)
add_executable(opt2 opt2.cpp)
add_executable(work work.cpp)
add_executable(scorer scorer.cpp)
add_executable(main main.cpp)
