
cmake_minimum_required(VERSION 3.2)
project(girgs)

set(CMAKE_CXX_STANDARD 14)
# set(CMAKE_VERBOSE_MAKEFILE on)

add_compile_options(
        #-Wall
        #-mavx2
)


SET(HEADER_FILES
        include/Generator.h
        include/SpatialTree.h
        include/SpatialTree.inl
        include/DebugHelper.h
        include/SpatialTreeCoordinateHelper.h
        include/SpatialTreeCoordinateHelper.inl
        include/WeightLayer.h
        include/WeightLayer.inl
        include/Common.h)


set(SOURCE_FILES
        source/Generator.cpp
        source/Common.cpp)


set(EXECUTABLE_FILES
        tools/benchmark.cpp
        tools/dev.cpp
        tools/test_avg_deg_estimation.cpp
        tools/test_Generator.cpp
        tools/test_SpartialTree.cpp
        )


include_directories(./include)
# add_library(girgs STATIC ${SOURCE_FILES} ${HEADER_FILES})
add_library(girgs SHARED ${SOURCE_FILES} ${HEADER_FILES})
foreach(file ${EXECUTABLE_FILES})
    get_filename_component(executable ${file} NAME_WE)
    add_executable(${executable} ${file})
    # target_link_libraries(${executable} girgs -static -static-libgcc -static-libstdc++)
    target_link_libraries(${executable} girgs)
endforeach()
