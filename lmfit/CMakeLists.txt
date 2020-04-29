include(GNUInstallDirs)

add_library(lmfit_local STATIC lmmin.c lmmin.h lmstruct.h lminvert.h)

# find_package(lmfit REQUIRED)
# target_link_libraries(centroids PUBLIC lmfit::lmfit)
# 
# find_package(OpenMP)
# if(OpenMP_FOUND)
# 	target_link_libraries(centroids PUBLIC OpenMP::OpenMP_CXX)
# endif()
# 
# if(CPPLINT_CHECK)
# 	add_dependencies(centroids ${CPPLINT_TARGET})
# 	message(STATUS "Checking CXX Code via cpplint")
# endif()
# 
# if (LIB_INSTALL)
# 	include(CMakePackageConfigHelpers)
# 
#     install(
#         TARGETS centroids 
# 		EXPORT centroids-targets
# 		LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
# 		PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
#         COMPONENT Libraries
# 	)
# 
# 	install(
# 		EXPORT centroids-targets
# 		FILE centroids-targets.cmake
# 		NAMESPACE centroids::
# 		DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/centroids"
# 	)
# 
# 	configure_package_config_file(
# 		"${CMAKE_CURRENT_SOURCE_DIR}/centroids-config.cmake.in"
# 		"${CMAKE_CURRENT_BINARY_DIR}/centroids-config.cmake"
# 		INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/centroids"
# 	)
# 	
# 	write_basic_package_version_file(
# 		"${CMAKE_CURRENT_BINARY_DIR}/centroids-config-version.cmake"
# 		VERSION ${centroids_VERSION}
# 		COMPATIBILITY AnyNewerVersion
# 	)
# 
# 	install(
# 		FILES "${CMAKE_CURRENT_BINARY_DIR}/centroids-config.cmake"
# 		"${CMAKE_CURRENT_BINARY_DIR}/centroids-config-version.cmake"
# 		DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/centroids"
# 	)
# endif()
