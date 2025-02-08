cmake_minimum_required(VERSION 3.5.1)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_PREFIX /usr/local)
project(tpsa)

#Bring the headers into the project
include_directories(include)

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.cc")

#Generate the shared library and the static library from the sources
add_library(tpsaShared SHARED ${SOURCES})
set_target_properties(tpsaShared PROPERTIES 
    OUTPUT_NAME "tpsa"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)
add_library(tpsaStatic STATIC ${SOURCES})
set_target_properties(tpsaStatic PROPERTIES 
    OUTPUT_NAME "tpsa"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

add_custom_target(build_examples
    COMMAND ${CMAKE_MAKE_PROGRAM} INSTALL_DIR=${CMAKE_INSTALL_PREFIX} all
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples
    COMMENT "Building examples using Makefile with custom INSTALL_DIR"
)

add_custom_target(build_tests
    COMMAND ${CMAKE_MAKE_PROGRAM} INSTALL_DIR=${CMAKE_INSTALL_PREFIX} all
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test
    COMMENT "Building tests using Makefile with custom INSTALL_DIR"
)

#Set the location for library installation -- i.e., /usr/local/lib in this case
#Use "sudo make install" to apply
install(TARGETS tpsaShared tpsaStatic DESTINATION)

#Install the header file
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/da.h DESTINATION include)

