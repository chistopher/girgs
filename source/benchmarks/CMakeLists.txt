
#
# Configure test project and environment
#

# Build google benchmark
download_project(PROJ                googlebenchmark
                 GIT_REPOSITORY      https://github.com/google/benchmark.git
                 GIT_TAG             master
                 UPDATE_DISCONNECTED 1
)

set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF)
add_subdirectory(
    ${googlebenchmark_SOURCE_DIR}
    ${googlebenchmark_BINARY_DIR})

# configure targets
foreach(target benchmark benchmark_main)
    set_target_properties(${target} PROPERTIES
        FOLDER "${IDE_FOLDER}"
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
        ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    )
endforeach()


# add onw benchmarks
add_subdirectory(bmi-benchmarks)
add_subdirectory(math-benchmarks)
