
#
# Configure test project and environment
#


# Get gtest
FetchContent_Declare(
  GTest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.14.0
  FIND_PACKAGE_ARGS
)

set(gtest_force_shared_crt  ON  CACHE BOOL "" FORCE) # Prevent GoogleTest from overriding our compiler/linker options when building with Visual Studio
set(INSTALL_GTEST OFF)
set(BUILD_GMOCK OFF) # see bug https://github.com/google/googletest/issues/4384
# since CMake 3.24 this tries a find_package(GTest) which should find GTest installed with `sudo apt install libgtest-dev`
# earlier CMake or a failed find_package(GTest) call lead to a download in the build dir
FetchContent_MakeAvailable(GTest) 

# windows needs gtest.dll in PATH. Installed versions should have this, but downloaded do not.
# To find DLL at test runtime (and during gtest_discover_tests()), we generate the gtest.dll in the same dir as the test .exe
# see https://github.com/google/googletest/issues/3565
if(WIN32)
  set_target_properties(gtest      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
  set_target_properties(gtest_main PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
endif()

add_subdirectory(girgs-test)
add_subdirectory(hypergirgs-test)
