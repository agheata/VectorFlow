#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
include_directories("${PROJECT_SOURCE_DIR}/event/inc")

file(GLOB sources src/*.cxx)
file(GLOB headers inc/*.h)

#----------------------------------------------------------------------------------------------
# Add the library
#
add_library(VFevent SHARED ${sources} )
target_link_libraries(VFevent -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECGEOM_LIBRARIES})

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
