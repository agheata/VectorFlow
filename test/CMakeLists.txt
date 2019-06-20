include_directories(
  ${CMAKE_SOURCE_DIR}/interface/inc
  ${CMAKE_SOURCE_DIR}/base/inc
)
#----------------------------------------------------------------------------------------------
# Add the executables
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

add_executable(testInterfaces testInterfaces.cxx)
target_link_libraries(testInterfaces TrackML -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECGEOM_LIBRARIES})
message(STATUS "Created target for testInterfaces")

#----------------------------------------------------------------------------------------------
# Install the files
#
install(TARGETS testInterfaces DESTINATION bin/tests)
