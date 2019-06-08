cmake_minimum_required(VERSION 3.12)
project(centroids)

set(CMAKE_CXX_STANDARD 14)

# Add a custom command that produces version.cpp, plus
# a dummy output that's not actually produced, in order
# to force version.cmake to always be re-run before the build
ADD_CUSTOM_COMMAND(
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
	${CMAKE_CURRENT_BINARY_DIR}/_version.cpp
	COMMAND ${CMAKE_COMMAND} -P
	${CMAKE_CURRENT_SOURCE_DIR}/version.cmake)

find_package(pybind11 REQUIRED)
pybind11_add_module(pycentroids 
	src/binding.cpp src/photons.cpp src/photons.h 
	${CMAKE_CURRENT_BINARY_DIR}/version.cpp
)
target_compile_options(pycentroids PRIVATE -Werror -pedantic -Wall -Wextra)

add_library(centroids 
	src/photons.cpp src/photons.h 
	${CMAKE_CURRENT_BINARY_DIR}/version.cpp
)
target_compile_options(centroids PRIVATE -Werror -pedantic -Wall -Wextra)
