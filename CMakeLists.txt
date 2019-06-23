cmake_minimum_required(VERSION 3.12)
project(centroids VERSION 0.1)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP REQUIRED)

include(PreventInSourceBuilds)

option(PEDANTIC "Compile with pedantic warnings" ON)
option(WERROR "Treat warnings as errors" OFF)
option(DEBUG_OUTPUT "Print debug outout to stderr" OFF)
option(LIB_INSTALL "Install library" ON)
option(BUILD_DOCS "Build documentation" ON)
option(CPPLINT_CHECK "Check CXX Source for Formatting" ON)
option(BUILD_PYTHON "Build Python Bindings using pybind11" ON)
option(BUILD_LIB "Build Library" ON)

if(CPPLINT_CHECK)
	include(cpplint)
	cpplint_add_subdirectory(lib)
	cpplint_add_subdirectory(python)
endif()

if(PEDANTIC)
	add_compile_options(-pedantic -Wall -Wextra)
endif()

if(WERROR)
	add_compile_options(-Werror)
endif()

if(DEBUG_OUTPUT)
 	add_compile_options(-DDEBUG_OUTPUT)
endif()

if(BUILD_PYTHON)
	add_subdirectory(python)
endif()

if(BUILD_DOCS)
	add_subdirectory(docs/sphinx)
	add_subdirectory(docs/doxygen)
endif()

if(BUILD_LIB)
	add_subdirectory(lib)
endif()

# Add __FILENAME__ with short path

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")

# This uses the lmfit lib

find_library(LMFIT_LIB lmfit)
find_path(LMFIT_INCLUDE lmmin.h)
message(STATUS "Found lmfit library: ${LMFIT_LIB}")
message(STATUS "Found lmmin.h : ${LMFIT_INCLUDE}")

# Add a custom command that produces version.cpp, plus
# a dummy output that's not actually produced, in order
# to force version.cmake to always be re-run before the build
add_custom_command(
	OUTPUT ${CMAKE_BINARY_DIR}/version.cpp
	${CMAKE_BINARY_DIR}/_version.cpp
	COMMAND ${CMAKE_COMMAND} -P
	${CMAKE_SOURCE_DIR}/cmake/version.cmake
)	
add_custom_target(version_info DEPENDS ${CMAKE_BINARY_DIR}/version.cpp)
