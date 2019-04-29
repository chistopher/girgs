
#
# Executable name and options
#

find_package(Boost 1.46)


# Target name
set(target hypergirgs-test)
message(STATUS "Test ${target}")


#
# Sources
#

set(sources
    main.cpp
    AngleHelper_test.cpp
    HyperbolicTree_test.cpp
    Point_test.cpp
    RadiusLayer_test.cpp
)


#
# Create executable
#

# Build executable
add_executable(${target}
    ${sources}
)

# Create namespaced alias
add_executable(${META_PROJECT_NAME}::${target} ALIAS ${target})


#
# Project options
#

set_target_properties(${target}
    PROPERTIES
    ${DEFAULT_PROJECT_OPTIONS}
    FOLDER "${IDE_FOLDER}"
)


#
# Include directories
#

target_include_directories(${target}
    PRIVATE
    ${DEFAULT_INCLUDE_DIRECTORIES}
    ${PROJECT_BINARY_DIR}/source/include
)

if(Boost_FOUND)
    message(STATUS "Found Boost. Including additional tests.")
    target_include_directories(${target} PRIVATE ${Boost_INCLUDE_DIRS})
    target_sources(${target} PRIVATE EdgeProbabilities_test.cpp)
else()
    message(STATUS "Boost not found. Skipping additional tests.")
endif()

#
# Libraries
#

target_link_libraries(${target}
    PRIVATE
    ${DEFAULT_LIBRARIES}
    ${META_PROJECT_NAME}::hypergirgs
    gmock
    gtest
)


#
# Compile definitions
#

target_compile_definitions(${target}
    PRIVATE
    ${DEFAULT_COMPILE_DEFINITIONS}
)


#
# Compile options
#

target_compile_options(${target}
    PRIVATE
    ${DEFAULT_COMPILE_OPTIONS}
)


#
# Linker options
#

target_link_libraries(${target}
    PRIVATE
    ${DEFAULT_LINKER_OPTIONS}
)
