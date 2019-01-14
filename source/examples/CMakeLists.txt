
# Check if examples are enabled
if(NOT OPTION_BUILD_EXAMPLES)
    return()
endif()

# Example applications
add_subdirectory(bench)
add_subdirectory(dev)
add_subdirectory(dev2)
add_subdirectory(dev3)
add_subdirectory(girggen)
add_subdirectory(hyper)
