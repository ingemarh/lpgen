cmake_minimum_required(VERSION 3.10)
set(CMAKE_C_COMPILER gcc)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_VERBOSE_MAKEFILE ON)
add_compile_options(-O4 -cpp -ffast-math -march=native)

project(dsp)

option(GUPTHREAD "Build with pthreads" ON)
#add_compile_options(-O4 -Wall -march=native)

add_library(plwin SHARED
)
target_sources(plwin PRIVATE
	plwin.c get_tid.c
)
target_link_libraries(plwin PUBLIC fftw3f pthread)


add_library(alt_decoder SHARED
)
target_sources(alt_decoder PRIVATE
	alt_decoder.c get_tid.c
)
target_link_libraries(alt_decoder PUBLIC pthread)


add_library(clutter SHARED
)
target_sources(clutter PRIVATE
	clutter.c get_tid.c
)
target_link_libraries(clutter PUBLIC fftw3f pthread)

add_library(pllow SHARED
)
target_sources(pllow PRIVATE
	pllow.c get_tid.c
)
target_link_libraries(pllow PUBLIC fftw3f pthread)

SET_TARGET_PROPERTIES(plwin clutter alt_decoder pllow PROPERTIES PREFIX "")
