cmake_minimum_required(VERSION 3.5.1)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_PREFIX /usr)
# Scikit-build-core sets these values for you, or you can just hard-code the
# name and version.
project(
  tpsa
  VERSION 1.0.0
  LANGUAGES CXX)
  
# Find the module development requirements (requires FindPython from 3.17 or
# scikit-build-core's built-in backport)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

#Set "TPSA_DIR" to the parent directory of the current source directory
get_filename_component(TPSA_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

#Bring the headers into the project
include_directories(${TPSA_DIR}/include)

#Add source files. The file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.cc" "${TPSA_DIR}/src/*.cc")

#message(STATUS "The value of TPSA_DIR is: ${TPSA_DIR}")
#message(STATUS "The value of PSOURCES is: ${SOURCES}")

# Add a library using FindPython's tooling (pybind11 also provides a helper like
# this)
python_add_library(_core MODULE ${SOURCES} WITH_SOABI)
target_link_libraries(_core PRIVATE pybind11::headers m)

# The install directory is the output (wheel) directory
install(TARGETS _core DESTINATION tpsa)
