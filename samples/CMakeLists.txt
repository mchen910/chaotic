project(chaotic LANGUAGES CXX)

set (CHAOTIC_SOURCE_FILES
	main.cpp)

set (CHAOTIC_HEADER_FILES
	../include/diffeq/diffeq.h
)

add_executable(samples ${CHAOTIC_SOURCE_FILES} ${CHAOTIC_HEADER_FILES})
target_include_directories(samples PRIVATE ../include/diffeq/)
target_link_libraries(samples PUBLIC chaotic ${OPENGL_LIBRARIES})

# target_link_libraries(samples PUBLIC chaotic glfw glad ${OPENGL_LIBRARIES})

