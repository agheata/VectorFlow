#----------------------------------------------------------------------------------------------
# Add the executables
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

include_directories(
  ${CMAKE_SOURCE_DIR}/base/inc
  ${CMAKE_SOURCE_DIR}/examples/trackML/inc
  ${CMAKE_SOURCE_DIR}/event/inc
  ${CMAKE_SOURCE_DIR}/magneticfield/inc
)

add_executable(testCocktailGenerator testCocktailGenerator.cxx)
target_link_libraries(testCocktailGenerator VFevent  ${VECGEOM_LIBRARIES})
message(STATUS "Created target for testCocktailGenerator")

add_executable(testTrackPropagation testTrackPropagation.cxx)
target_link_libraries(testTrackPropagation VFevent Vmagfield TrackML  ${VECGEOM_LIBRARIES})
message(STATUS "Created target for testTrackPropagation")

#----------------------------------------------------------------------------------------------
# Install the files
#
install(TARGETS testCocktailGenerator testTrackPropagation DESTINATION bin/tests)
