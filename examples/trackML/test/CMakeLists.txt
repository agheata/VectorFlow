include_directories(${CMAKE_SOURCE_DIR}/examples/trackML/inc)
#----------------------------------------------------------------------------------------------
# Add the executables
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

add_executable(testTrackMLgeometry testTrackMLgeometry.cxx)
target_link_libraries(testTrackMLgeometry TrackML VFevent -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECGEOM_LIBRARIES})
message(STATUS "Created target for testTrackMLgeometry")

#----------------------------------------------------------------------------------------------
# Install the files
#
install(TARGETS testTrackMLgeometry DESTINATION bin/tests)
