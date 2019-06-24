include_directories(
  ${CMAKE_SOURCE_DIR}/interface/inc
  ${CMAKE_SOURCE_DIR}/base/inc
  ${CMAKE_SOURCE_DIR}/event/inc
  ${CMAKE_SOURCE_DIR}/magneticfield/inc
  ${CMAKE_SOURCE_DIR}/examples/trackML/inc
)
#----------------------------------------------------------------------------------------------
# Add the executables
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

add_executable(testInterfaces testInterfaces.cxx)
target_link_libraries(testInterfaces TrackML -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECGEOM_LIBRARIES})
message(STATUS "Created target for testInterfaces")

add_executable(testTrackPropagationPipeLine testTrackPropagationPipeLine.cxx)
target_link_libraries(testTrackPropagationPipeLine TrackML VFevent Vmagfield -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECGEOM_LIBRARIES})
message(STATUS "Created target for testTrackPropagationPipeLine")

#----------------------------------------------------------------------------------------------
# Install the files
#
install(TARGETS testInterfaces DESTINATION bin/tests)

install(TARGETS testTrackPropagationPipeLine DESTINATION bin/tests)
