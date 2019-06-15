cmake_minimum_required(VERSION 3.12)
project(centroids VERSION 0.1)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(OpenMP REQUIRED)

include(PreventInSourceBuilds)

include(cpplint)
cpplint_add_subdirectory(lib)
cpplint_add_subdirectory(python)

set(CMAKE_CXX_STANDARD 14)

option(PEDANTIC "Compile with pedantic warnings" ON)
option(WERROR "Treat warnings as errors" OFF)
option(DEBUG_OUTPUT "Print debug outout to stderr" OFF)
option(LIB_INSTALL "Install library" ON)

if(PEDANTIC)
	add_compile_options(-pedantic -Wall -Wextra)
endif()

if(WERROR)
	add_compile_options(-Werror)
endif()

if(DEBUG_OUTPUT)
 	add_compile_options(-DDEBUG_OUTPUT)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")

# Add a custom command that produces version.cpp, plus
# a dummy output that's not actually produced, in order
# to force version.cmake to always be re-run before the build
add_custom_command(
	OUTPUT ${CMAKE_BINARY_DIR}/version.cpp
	${CMAKE_BINARY_DIR}/_version.cpp
	COMMAND ${CMAKE_COMMAND} -P
	${CMAKE_SOURCE_DIR}/cmake/version.cmake
	DEPENDS ${CMAKE_BINARY_DIR}/_version.cpp
)	

add_custom_target(version_info DEPENDS ${CMAKE_BINARY_DIR}/version.cpp)

add_subdirectory(lib)
add_subdirectory(python)

