
#
# Configure test project and environment
#

# Build google benchmark
FetchContent_Declare(
  benchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG        v1.8.3
  FIND_PACKAGE_ARGS
)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(benchmark) 

# add onw benchmarks
add_subdirectory(bmi-benchmarks)
add_subdirectory(math-benchmarks)
