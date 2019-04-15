
# Check if examples are enabled
if(NOT OPTION_BUILD_EXAMPLES)
    return()
endif()

# Example applications
add_subdirectory(bench)
add_subdirectory(bench_hyper)
add_subdirectory(dev3)
add_subdirectory(dev-girgs)
add_subdirectory(girggen)
add_subdirectory(hyper)
add_subdirectory(hyper_range)
