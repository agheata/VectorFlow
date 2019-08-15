include_directories(
  ${CMAKE_SOURCE_DIR}/interface/inc
  ${CMAKE_SOURCE_DIR}/base/inc
  ${CMAKE_SOURCE_DIR}/event/inc
  ${CMAKE_SOURCE_DIR}/magneticfield/inc
  ${CMAKE_SOURCE_DIR}/examples/trackML/inc
  ${CMAKE_SOURCE_DIR}/test/gbench/include
)
#----------------------------------------------------------------------------------------------
# Add the executables
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

add_executable(testInterfaces testInterfaces.cxx)
target_link_libraries(testInterfaces TrackML ${VECGEOM_LIBRARIES})
message(STATUS "Created target for testInterfaces")

add_executable(testTrackPropagationPipeLine testTrackPropagationPipeLine.cxx)
target_link_libraries(testTrackPropagationPipeLine TrackML VFevent Vmagfield ${VECGEOM_LIBRARIES})
message(STATUS "Created target for testTrackPropagationPipeLine")

add_executable(testVectorTypesTrackPropagationPipeLine testVectorTypesTrackPropagationPipeLine.cxx)
target_link_libraries(testVectorTypesTrackPropagationPipeLine TrackML VFevent Vmagfield benchmark ${VECGEOM_LIBRARIES} ${PROFILE_LIBRARIES})
message(STATUS "Created target for testVectorTypesTrackPropagationPipeLine")

add_executable(testTrackPropagationComplex testTrackPropagationComplex.cxx)
target_link_libraries(testTrackPropagationComplex TrackML VFevent Vmagfield ${VECGEOM_LIBRARIES} ${PROFILE_LIBRARIES})
message(STATUS "Created target for testTrackPropagationComplex")

#----------------------------------------------------------------------------------------------
# Install the files
#
install(TARGETS testInterfaces DESTINATION bin/tests)
install(TARGETS testTrackPropagationPipeLine DESTINATION bin/tests)
install(TARGETS testVectorTypesTrackPropagationPipeLine DESTINATION bin/tests)
install(TARGETS testTrackPropagationComplex DESTINATION bin/tests)

add_subdirectory(gbench EXCLUDE_FROM_ALL)
