cmake_minimum_required(VERSION 3.6.0)

include(cmake/RecordCmdLine.cmake)

option(BUILD_TESTING  "Build test binaries and create test target"  ON)
option(BUILD_EXAMPLES "Build binaries for performance benchmarking" ON)
option(ROOT "Include ROOT" OFF)
set(VecCore_VERSION 0.5.0)

if(BUILD_TESTING)
  message("Vc backend selected for testing")
  set(VecCore_BACKEND Vc)
endif()

project(VectorFlow VERSION 0.1.0)

add_library(VectorFlow INTERFACE)
target_include_directories(VectorFlow INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)

# Added due to VecGeom dependency with ROOT
if ( ROOT )
  find_package(ROOT)
  # Added the following IF statement to NOT force ROOT package to be required
  if ( ROOT_FOUND )
    message(STATUS "ROOT dir found: ${ROOT_INCLUDE_DIRS}")
    message(STATUS "ROOT libraries: ${ROOT_LIBRARIES}")
    include_directories(${ROOT_INCLUDE_DIRS})
    target_link_libraries(VectorFlow INTERFACE ${ROOT_LIBRARIES})
  endif ( ROOT_FOUND )
endif ( ROOT )

#Find VecCore with correct backend
find_package(VecCore ${VecCore_VERSION} REQUIRED COMPONENTS ${VecCore_BACKEND})
target_link_libraries(VectorFlow INTERFACE VecCore::VecCore)

# Find VecMath library with correct backend
find_package(VecMath 0.1.0 REQUIRED COMPONENTS Vc)
#Temporary solution is to use VecMath_INCLUDES/LIBRARIES variables since GV does not support modern cmake targets
include_directories(${VecMath_INCLUDE_DIRS})
add_definitions(${VecMath_DEFINITIONS})

# Find VecGeom geometry headers library
find_package(VecGeom REQUIRED)
message(STATUS "VecGeom libraries ${VECGEOM_LIBRARIES}")
#---ensure that VectorFlow is compiled with same architecture flags as VecGeom
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VECGEOM_CXX_FLAGS}")
include_directories(${VECGEOM_INCLUDE_DIRS})
set ( VECGEOM_ROOT OFF )
if ( VECGEOM_DEFINITIONS MATCHES VECGEOM_ROOT )
  set ( VECGEOM_ROOT ON )
endif()

#installation configuration
install(DIRECTORY include/ DESTINATION include)

include(CMakePackageConfigHelpers)

#exported targets file
install(TARGETS VectorFlow EXPORT VectorFlowTargets
  INCLUDES DESTINATION include
)
export(EXPORT VectorFlowTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/VectorFlowTargets.cmake"
  NAMESPACE VectorFlow::
)

#package config file
set(VectorFlow_INSTALL_INCLUDEDIR include)
configure_package_config_file(cmake/VectorFlowConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/VectorFlowConfig.cmake
  INSTALL_DESTINATION lib/cmake/VectorFlow
  PATH_VARS VectorFlow_INSTALL_INCLUDEDIR
)

#version file
write_basic_package_version_file("VectorFlowConfigVersion.cmake"
	VERSION ${VectorFlow_VERSION}
	COMPATIBILITY SameMajorVersion
	)

#install package
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/VectorFlowConfigVersion.cmake
${CMAKE_CURRENT_BINARY_DIR}/VectorFlowConfig.cmake
${CMAKE_CURRENT_BINARY_DIR}/VectorFlowTargets.cmake
	DESTINATION lib/cmake/VectorFlow
	)

#----------------------------------------------------------------------------------------------
# Use absolute RPATH when linking to avoid needing to export LD_LIBRARY_PATH
#
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

add_subdirectory(base)
add_subdirectory(interface)

if (BUILD_TESTING)
  add_subdirectory(test)
  add_subdirectory(magneticfield)
endif()


if (BUILD_EXAMPLES)
  add_subdirectory(event)
  add_subdirectory(examples)
endif()
