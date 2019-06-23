set(lib centroids)
set(${lib}_LIBRARY ${lib} PARENT_SCOPE)

set_property(SOURCE ${CMAKE_BINARY_DIR}/version.cpp PROPERTY GENERATED 1)
set(src_files photons.cpp ${CMAKE_BINARY_DIR}/version.cpp)
set(inc_files centroids.h)

add_library(${lib} SHARED ${src_files})

add_dependencies(${lib} version_info)

target_link_libraries(${lib} PUBLIC "${LMFIT_LIB}")
target_link_libraries(${lib} PUBLIC OpenMP::OpenMP_CXX)
target_include_directories(${lib} PUBLIC ${LMFIT_INCLUDE})

set_target_properties(
    ${lib} PROPERTIES
    OUTPUT_NAME ${lib}
    VERSION ${centroids_VERSION_MAJOR}
    SOVERSION ${centroids_VERSION_MINOR})

if(CPPLINT_CHECK)
	add_dependencies(${lib} ${CPPLINT_TARGET})
	message(STATUS "Checking CXX Code via cpplint")
endif()

if (LIB_INSTALL)
    install(
        TARGETS ${lib} LIBRARY
		DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        COMPONENT Libraries)
    install(
        FILES ${inc_files}
		DESTINATION ${CMAKE_INSTALL_PREFIX}/include
        COMPONENT Headers)
endif()
