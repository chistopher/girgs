
# Check if CLI's are enabled
if(NOT OPTION_BUILD_CLI)
    return()
endif()

# Example applications
add_subdirectory(gengirg)
add_subdirectory(genhrg)
