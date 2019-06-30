include(GNUInstallDirs)
set(centroids_LIBRARY centroids PARENT_SCOPE)
set(SOURCE_FILES photons.cpp ${CMAKE_BINARY_DIR}/version.cpp)
set(PUBLIC_INCLUDE_FILES centroids.h)
set(PRIVATE_INCLUDE_FILES photons.h)

set_property(SOURCE ${CMAKE_BINARY_DIR}/version.cpp PROPERTY GENERATED 1)

add_library(centroids SHARED ${SOURCE_FILES} 
	${PUBLIC_INCLUDE_FILES} ${PRIVATE_INCLUDE_FILES})

set_target_properties(
    centroids PROPERTIES
    OUTPUT_NAME centroids
    VERSION ${centroids_VERSION}
	SOVERSION ${centroids_VERSION_MAJOR}
	PUBLIC_HEADER "${PUBLIC_HEADER_FILES}"
	PRIVATE_HEADER "${PRIVATE_HEADER_FILES}"
)

target_include_directories(centroids PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

add_dependencies(centroids version_info)

find_package(lmfit REQUIRED)
target_link_libraries(centroids PUBLIC lmfit::lmfit)

find_package(OpenMP)
if(OpenMP_FOUND)
	target_link_libraries(centroids PUBLIC OpenMP::OpenMP_CXX)
endif()

if(CPPLINT_CHECK)
	add_dependencies(centroids ${CPPLINT_TARGET})
	message(STATUS "Checking CXX Code via cpplint")
endif()

if (LIB_INSTALL)
	include(CMakePackageConfigHelpers)

    install(
        TARGETS centroids 
		EXPORT centroids-targets
		LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
		PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        COMPONENT Libraries
	)

	install(
		EXPORT centroids-targets
		FILE centroids-targets.cmake
		NAMESPACE centroids::
		DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/centroids"
	)

	configure_package_config_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/centroids-config.cmake.in"
		"${CMAKE_CURRENT_BINARY_DIR}/centroids-config.cmake"
		INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/centroids"
	)
	
	write_basic_package_version_file(
		"${CMAKE_CURRENT_BINARY_DIR}/centroids-config-version.cmake"
		VERSION ${centroids_VERSION}
		COMPATIBILITY AnyNewerVersion
	)

	install(
		FILES "${CMAKE_CURRENT_BINARY_DIR}/centroids-config.cmake"
		"${CMAKE_CURRENT_BINARY_DIR}/centroids-config-version.cmake"
		DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/centroids"
	)
endif()
