cmake_minimum_required(VERSION 3.8)

project(chaotic LANGUAGES CXX)

add_subdirectory(src)

option(CHAOTIC_BUILD_SAMPLES "Build the sample program" ON)

if (CHAOTIC_BUILD_SAMPLES)

	add_subdirectory(extern/glad)
	add_subdirectory(extern/glfw)
	add_subdirectory(samples)

	# default startup project for Visual Studio
	if (MSVC)
		set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT samples)
	endif()

endif()
