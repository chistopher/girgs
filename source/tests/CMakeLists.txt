
#
# Configure test project and environment
#


# Build gtest/gmock
download_project(PROJ                googletest
                 GIT_REPOSITORY      https://github.com/google/googletest.git
                 GIT_TAG             main
                 UPDATE_DISCONNECTED 1
)

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt  ON  CACHE BOOL "" FORCE)
set(INSTALL_GTEST           OFF CACHE BOOL "")
set(gtest_build_samples     OFF CACHE BOOL "")
set(gtest_build_tests       OFF CACHE BOOL "")
set(gmock_build_tests       OFF CACHE BOOL "")

add_subdirectory(
    ${googletest_SOURCE_DIR}
    ${googletest_BINARY_DIR}
)

# configure targets
foreach(target gtest gtest_main gmock gmock_main)
    set_target_properties(${target} PROPERTIES
        FOLDER "${IDE_FOLDER}"
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
        ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
        PDB_OUTPUT_DIRECTORY     ${PROJECT_BINARY_DIR}
        CXX_STANDARD 11 # somehow needed in ubuntu 16.04
    )
endforeach()

add_subdirectory(girgs-test)
add_subdirectory(hypergirgs-test)
