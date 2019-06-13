cmake_minimum_required(VERSION 3.12)
project(centroids)

include(cmake/cpplint.cmake)
cpplint_add_subdirectory(src)

set(CMAKE_CXX_STANDARD 14)
option(PEDANTIC "Compile with pedantic warnings" ON)
option(WERROR "Treat warnings as errors" OFF)
option(DEBUG_OUTPUT "Print debug outout to stderr" ON)

if(PEDANTIC)
	add_compile_options(-pedantic -Wall -Wextra)
endif()

if(WERROR)
	add_compile_options(-Werror)
endif()

if(DEBUG_OUTPUT)
	add_compile_options("-DDEBUG_OUTPUT")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")

# Add a custom command that produces version.cpp, plus
# a dummy output that's not actually produced, in order
# to force version.cmake to always be re-run before the build
ADD_CUSTOM_COMMAND(
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
	${CMAKE_CURRENT_BINARY_DIR}/_version.cpp
	COMMAND ${CMAKE_COMMAND} -P
	${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.cmake)

find_package(pybind11 REQUIRED)
pybind11_add_module(pycentroids
	src/binding.cpp src/photons.cpp src/photons.h
	${CMAKE_CURRENT_BINARY_DIR}/version.cpp
)

add_library(centroids
	src/photons.cpp src/photons.h
	${CMAKE_CURRENT_BINARY_DIR}/version.cpp
)
