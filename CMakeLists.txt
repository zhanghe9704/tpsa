cmake_minimum_required(VERSION 3.5.1)
project(tpsa)
set(CMAKE_BUILD_TYPE Release)
set (CMAKE_CXX_STANDARD 14)

#Bring the headers, such as Student.h into the project
include_directories(include)

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.cc")

#Generate the shared library from the sources
add_library(tpsaso SHARED ${SOURCES})
add_library(tpsa STATIC ${SOURCES})

#Set the location for library installation -- i.e., /usr/lib in this case
# not really necessary in this example. Use "sudo make install" to apply
#install(TARGETS libtpsa DESTINATION /usr/lib)
