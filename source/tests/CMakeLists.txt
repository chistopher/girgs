
#
# Configure test project and environment
#


# Get gtest/gmock
include(FetchContent)
FetchContent_Declare(
  GTest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.14.0
  FIND_PACKAGE_ARGS
)


# # Prevent GoogleTest from overriding our compiler/linker options when building with Visual Studio
set(gtest_force_shared_crt  ON  CACHE BOOL "" FORCE)
# since CMake 3.24 this tries a find_package(GTest) which should find GTest installed with `sudo apt install libgtest-dev`
# earlier CMake or a failed find_package(GTest) call lead to a download in the build dir
FetchContent_MakeAvailable(GTest) 

add_subdirectory(girgs-test)
add_subdirectory(hypergirgs-test)
