include_directories(
  ${CMAKE_SOURCE_DIR}/examples/trackML/inc
  ${CMAKE_SOURCE_DIR}/magneticfield/inc
  ${CMAKE_SOURCE_DIR}/event/inc
)
#----------------------------------------------------------------------------------------------
# Add the executables
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

add_executable(testTrackMLgeometry testTrackMLgeometry.cxx)
target_link_libraries(testTrackMLgeometry TrackML VFevent -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECGEOM_LIBRARIES})
message(STATUS "Created target for testTrackMLgeometry")

add_executable(testFieldManager testFieldManager.cxx)
target_link_libraries(testFieldManager TrackML VFevent Vmagfield -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECGEOM_LIBRARIES})
message(STATUS "Created target for testFieldManager")

#----------------------------------------------------------------------------------------------
# Install the files
#
install(TARGETS testTrackMLgeometry testFieldManager DESTINATION bin/tests)
