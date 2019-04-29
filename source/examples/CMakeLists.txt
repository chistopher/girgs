
# Check if examples are enabled
if(NOT OPTION_BUILD_EXAMPLES)
    return()
endif()

# Example applications
add_subdirectory(bench_girg)
add_subdirectory(bench_hyper)
add_subdirectory(inclusion_single)
add_subdirectory(inclusion_range)
