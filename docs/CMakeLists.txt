
if(NOT OPTION_BUILD_DOCS)
    return()
endif()

#
# Find doxygen
#

find_package(Doxygen)
if(NOT DOXYGEN_FOUND)
    message(STATUS "Disabled generation of doxygen documentation (missing doxygen).")
    return()
endif()

message(STATUS "Doc ${target}")


#
# Customize Doxygen
#

set(DOXYGEN_QUIET                   YES)
set(DOXYGEN_USE_MDFILE_AS_MAINPAGE  ${PROJECT_SOURCE_DIR}/README.md)
set(DOXYGEN_PROJECT_LOGO            ${PROJECT_SOURCE_DIR}/girgs-logo.png)
set(DOXYGEN_USE_MATHJAX             YES)
set(DOXYGEN_EXTRACT_ALL             YES)

#
# Create target
#

doxygen_add_docs(docs
    ${PROJECT_SOURCE_DIR}/README.md
    ${PROJECT_SOURCE_DIR}/source/girgs/include/girgs/Generator.h
    ${PROJECT_SOURCE_DIR}/source/hypergirgs/include/hypergirgs/Generator.h
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/source
)

#
# Deployment
#

install(
    DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
    DESTINATION ${INSTALL_DOC}
    COMPONENT docs
)
