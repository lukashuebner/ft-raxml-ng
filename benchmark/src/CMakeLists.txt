
message (STATUS "Building benchmarks")

file (GLOB_RECURSE RAXML_BENCHMARK_SOURCES ${PROJECT_SOURCE_DIR}/benchmark/src/*.cpp ${RAXML_SOURCES})

# sources list now has 2 Main.cpp, old has to be removed
list(REMOVE_ITEM RAXML_BENCHMARK_SOURCES "${PROJECT_SOURCE_DIR}/src/main.cpp")

include_directories (${PROJECT_SOURCE_DIR})

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/benchmark/bin)

add_executable        (raxml_benchmark_module ${RAXML_BENCHMARK_SOURCES})

target_link_libraries (raxml_benchmark_module ${RAXML_LIBS})

#if(GMP_FOUND)
#  target_link_libraries(raxml_benchmark_module ${GMP_LIBRARIES})
#endif()

#target_link_libraries (raxml_benchmark_module ${MPI_CXX_LIBRARIES})

#if(MPI_COMPILE_FLAGS)
#  set_target_properties(raxml_benchmark_module PROPERTIES
#  COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
#endif()

#if(MPI_LINK_FLAGS)
#  set_target_properties(raxml_benchmark_module PROPERTIES
#    LINK_FLAGS "${MPI_LINK_FLAGS}")
#endif()

set_target_properties (raxml_benchmark_module PROPERTIES OUTPUT_NAME raxml_benchmark)
set_target_properties (raxml_benchmark_module PROPERTIES PREFIX "")


#add_test (NAME raxml_test
#          COMMAND ${PROJECT_SOURCE_DIR}/test/bin/raxml_test
#          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
