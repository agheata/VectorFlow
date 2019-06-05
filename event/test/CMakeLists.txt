#----------------------------------------------------------------------------------------------
# Add the executables
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

add_executable(testCocktailGenerator testCocktailGenerator.cxx)
target_link_libraries(testCocktailGenerator VFevent -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECGEOM_LIBRARIES})
message(STATUS "Created target for testCocktailGenerator")

#----------------------------------------------------------------------------------------------
# Install the files
#
install(TARGETS testCocktailGenerator DESTINATION bin/tests)
