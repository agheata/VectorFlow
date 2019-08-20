#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
include_directories(
  "${PROJECT_SOURCE_DIR}/event/inc"
  "${PROJECT_SOURCE_DIR}/base/inc"
)

file(GLOB sources src/*.cxx)
file(GLOB headers inc/*.h)

#----------------------------------------------------------------------------------------------
# Add the library
#
add_library(VFevent SHARED ${sources} )
target_link_libraries(VFevent  ${VECGEOM_LIBRARIES} ${VecCore_LIBRARIES})

#----------------------------------------------------------------------------------------------
# Install the files
#
install(TARGETS VFevent LIBRARY DESTINATION lib)
install(FILES ${headers} DESTINATION inc)

if (BUILD_TESTING)
  add_subdirectory(test)
endif()

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
