find_package(pybind11 REQUIRED)

set_property(SOURCE ${CMAKE_BINARY_DIR}/version.cpp PROPERTY GENERATED 1)

pybind11_add_module(pycentroids EXCLUDE_FROM_ALL
	binding.cpp
	${CMAKE_SOURCE_DIR}/lib/photons.cpp
	${CMAKE_BINARY_DIR}/version.cpp
)

target_include_directories(pycentroids PUBLIC ${CMAKE_SOURCE_DIR}/lib)
add_dependencies(pycentroids version_info)

target_link_libraries(pycentroids PUBLIC liblmfit.so)
target_link_libraries(pycentroids PUBLIC OpenMP::OpenMP_CXX)
