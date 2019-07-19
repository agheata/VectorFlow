#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources src/*.cxx)
file(GLOB headers inc/*.h)
file(GLOB data data/*.csv)

include_directories(
  inc
  ${CMAKE_SOURCE_DIR}/base/inc
  ${CMAKE_SOURCE_DIR}/examples/trackML/inc
  ${CMAKE_SOURCE_DIR}/event/inc
  ${CMAKE_SOURCE_DIR}/magneticfield/inc
)

#----------------------------------------------------------------------------------------------
# Add the library
#
add_library(TrackML SHARED ${sources} )
target_link_libraries(TrackML VFevent Vmagfield  ${VECGEOM_LIBRARIES})

#----------------------------------------------------------------------------------------------
# Install the files
#
install(TARGETS TrackML LIBRARY DESTINATION lib)
install(FILES ${headers} DESTINATION examples/TrackML/inc)
install(FILES ${data} DESTINATION bin/tests/data)

add_subdirectory(test)

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()
