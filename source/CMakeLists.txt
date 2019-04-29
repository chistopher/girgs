
#
# Configuration for all sub-projects
#

# Generate version-header
configure_file(version.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/${META_PROJECT_NAME}/${META_PROJECT_NAME}-version.h)


#
# Sub-projects
#

# Libraries
set(IDE_FOLDER "")
add_subdirectory(girgs)
add_subdirectory(hypergirgs)

# Examples
set(IDE_FOLDER "Examples")
add_subdirectory(examples)

# Generator CLI's
set(IDE_FOLDER "CLI")
add_subdirectory(cli)

# Tests
if(OPTION_BUILD_TESTS AND NOT MINGW)
    set(IDE_FOLDER "Tests")
    add_subdirectory(tests)
endif()

# Benchmarks
if(OPTION_BUILD_BENCHMARKS)
    set(IDE_FOLDER "Benchmarks")
    add_subdirectory(benchmarks)
endif()


#
# Deployment
#

# Deploy generated headers
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/${META_PROJECT_NAME} DESTINATION include COMPONENT dev)
