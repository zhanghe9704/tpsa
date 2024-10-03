# Set the minimum CMake version
cmake_minimum_required(VERSION 3.5.1)

# Set the project name
project(TPSAJuliaWrapper)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Use Julia to find the CxxWrap prefix path
message(STATUS "Finding CxxWrap prefix path using Julia...")

# Run Julia to get CxxWrap.prefix_path() and capture the output
execute_process(
    COMMAND julia -E "using CxxWrap; print(CxxWrap.prefix_path())"
    OUTPUT_VARIABLE CXXWRAP_PREFIX_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Manually trim "nothing" if it's present
string(REPLACE "nothing" "" CXXWRAP_PREFIX_PATH "${CXXWRAP_PREFIX_PATH}")
string(STRIP CXXWRAP_PREFIX_PATH "${CXXWRAP_PREFIX_PATH}")

message(STATUS "Found CxxWrap prefix path: ${CXXWRAP_PREFIX_PATH}")

# Set the CMAKE_PREFIX_PATH to the CxxWrap prefix path
set(CMAKE_PREFIX_PATH "${CXXWRAP_PREFIX_PATH}/lib/cmake/JlCxx")

# Find CxxWrap (JlCxx)
find_package(JlCxx REQUIRED)

# Find CxxWrap, which provides the necessary CMake files
set(CMAKE_PREFIX_PATH "$ENV{HOME}/.julia/packages/CxxWrap/...") # Adjust this path

message(STATUS "Julia_INCLUDE_DIR: ${Julia_INCLUDE_DIR}")
message(STATUS "Julia_LIBRARY: ${Julia_LIBRARY}")
message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
message(STATUS "CxxWrap_INCLUDE_DIRS: ${CxxWrap_INCLUDE_DIRS}")

# Include the CxxWrap and Julia include directories
include_directories(${CMAKE_SOURCE_DIR}/../include)
include_directories(${JULIA_INCLUDE_DIR})
include_directories(${CxxWrap_INCLUDE_DIRS})

# Set the source files
file(GLOB SOURCES "{CMAKE_SOURCE_DIR}/../src/*.cc")
list(APPEND SOURCES
    ${CMAKE_SOURCE_DIR}/julia-wrapper.cc
)

# Add the library to build the shared object (libtpsajuliawrapper.so)
add_library(juliaTPSA SHARED ${SOURCES})

# Link the Julia library
target_link_libraries(juliaTPSA ${Julia_LIBRARY}  JlCxx::cxxwrap_julia)

# Install the shared library
install(TARGETS juliaTPSA DESTINATION lib)

