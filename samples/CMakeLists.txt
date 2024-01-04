cmake_minimum_required(VERSION 3.20)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories(${CMAKE_HOME_DIRECTORY}/include/)
add_executable(main main.cpp ${CMAKE_HOME_DIRECTORY}/include/diffeq.h)
target_compile_features(main PRIVATE cxx_std_17)

target_link_libraries(main PRIVATE glad glfw)