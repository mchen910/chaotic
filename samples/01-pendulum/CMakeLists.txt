cmake_minimum_required(VERSION 3.20)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(HEADER_FILES ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/samples)

add_executable(01-pendulum pendulum.cpp)
target_include_directories(01-pendulum PRIVATE ${HEADER_FILES})
target_compile_features(01-pendulum PRIVATE cxx_std_17)
target_link_libraries(01-pendulum PRIVATE glad glm glfw)