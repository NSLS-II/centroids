find_package(pybind11 REQUIRED)

set_property(SOURCE ${CMAKE_BINARY_DIR}/version.cpp PROPERTY GENERATED 1)

pybind11_add_module(_pycentroids EXCLUDE_FROM_ALL
	binding.cpp
	${CMAKE_SOURCE_DIR}/lib/photons.cpp
	${CMAKE_BINARY_DIR}/version.cpp
)

target_include_directories(_pycentroids PUBLIC
	${CMAKE_SOURCE_DIR}/lib)
target_include_directories(_pycentroids PUBLIC
	${LMFIT_INCLUDE})

add_dependencies(_pycentroids version_info)

if(CPPLINT_CHECK)
	add_dependencies(_pycentroids ${CPPLINT_TARGET})
	message(STATUS "Checking CXX Code via cpplint")
endif()

target_link_libraries(_pycentroids PUBLIC ${LMFIT_LIB})
if(OpenMP_FOUND)
	target_link_libraries(_pycentroids PUBLIC OpenMP::OpenMP_CXX)
endif()