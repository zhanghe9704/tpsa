cmake_minimum_required(VERSION 3.5.1)

set(CMAKE_C_COMPILER "/usr/bin/clang")
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_PREFIX /usr)
project(tpsa)

#Bring the headers into the project
include_directories(include)

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.cc")

#Generate the shared library and the static llibary from the sources
add_library(tpsaShared SHARED ${SOURCES})
set_target_properties(tpsaShared PROPERTIES OUTPUT_NAME "tpsa")
add_library(tpsaStatic STATIC ${SOURCES})
set_target_properties(tpsaStatic PROPERTIES OUTPUT_NAME "tpsa")

#Set the location for library installation -- i.e., /usr/lib in this case
#Use "sudo make install" to apply
install(TARGETS tpsaShared tpsaStatic DESTINATION)

#Install the header file
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/da.h DESTINATION include)

