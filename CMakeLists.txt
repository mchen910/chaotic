cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_COMPILER "/opt/homebrew/Cellar/gcc/11.2.0_3/bin/gcc-11")
project(chaotic LANGUAGES CXX)

add_subdirectory(extern/glad)
add_subdirectory(extern/glfw)
add_subdirectory(samples)
