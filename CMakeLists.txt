cmake_minimum_required(VERSION 3.10)
project(DNA_encoding)
include_directories(./include)
aux_source_directory(./src DIR_SRCS)
set(CMAKE_CXX_STANDARD 11)

add_executable(DNA_encoding main.cpp src/global.cpp src/DNA_encoder.cpp include/DNA_encoder.h src/rabin.cpp include/rabin.h)