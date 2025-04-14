cmake_minimum_required(VERSION 3.22)
project(ePIC++ LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -O3")

# Source files for main executable
file(GLOB SOURCES "src/*.cpp")

# Find Python with Development and NumPy components
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
if(NOT Python3_FOUND)
    message(FATAL_ERROR "Python3/Python3-dev not found. Please install Python3/Python3-dev")
endif()

# Find HDF5 with CXX component
find_package(HDF5 COMPONENTS CXX REQUIRED)
if(NOT HDF5_FOUND)
    message(FATAL_ERROR "HDF5 not found. Please install HDF5.")
endif()

# Find Threads
find_package(Threads REQUIRED)

# Find Eigen (optional)
find_package(Eigen3 QUIET)
if(EIGEN3_FOUND)
    add_definitions(-DEIGEN_AVAILABLE)
    include_directories(${EIGEN3_INCLUDE_DIR})
endif()

# Print version information
message(STATUS "Found Python version: ${Python3_VERSION}")
message(STATUS "Python executable: ${Python3_EXECUTABLE}")
message(STATUS "Python include dirs: ${Python3_INCLUDE_DIRS}")
message(STATUS "NumPy include dirs: ${Python3_NumPy_INCLUDE_DIRS}")
message(STATUS "Found HDF5 version: ${HDF5_VERSION}")
message(STATUS "HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")
if(EIGEN3_FOUND)
    message(STATUS "Found Eigen3 at: ${EIGEN3_INCLUDE_DIR}")
else()
    message(STATUS "Eigen3 not found; Eigen solvers will be skipped.")
endif()

# Include directories
include_directories(
    include
    linearalgebra
    ${HDF5_INCLUDE_DIRS}
    ${Python3_INCLUDE_DIRS}
    ${Python3_NumPy_INCLUDE_DIRS}
)

# Main executable
add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} PRIVATE
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    Threads::Threads
    ${Python3_LIBRARIES}
    Python3::NumPy  # Use CMake-provided target for NumPy
    m
)
target_include_directories(${PROJECT_NAME} PRIVATE
    ${HDF5_INCLUDE_DIRS}
    ${Python3_INCLUDE_DIRS}
    ${Python3_NumPy_INCLUDE_DIRS}
)

# Enable testing
enable_testing()

# Test executable for solver performance
add_executable(test_solvers tests/test_solvers.cpp)
target_link_libraries(test_solvers PRIVATE
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    Threads::Threads
    ${Python3_LIBRARIES}
    Python3::NumPy
    m
)
if(EIGEN3_FOUND)
    target_link_libraries(test_solvers PRIVATE Eigen3::Eigen)
endif()
target_include_directories(test_solvers PRIVATE
    include
    linearalgebra
    ${HDF5_INCLUDE_DIRS}
    ${Python3_INCLUDE_DIRS}
    ${Python3_NumPy_INCLUDE_DIRS}
)

# Add test to CTest
add_test(NAME SolverPerformanceTest COMMAND test_solvers)

# Custom target to run tests
add_custom_target(run_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    DEPENDS test_solvers
    COMMENT "Running solver performance tests"
)


add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/CMakeCache.txt
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/cmake_install.cmake
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/Makefile
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/${PROJECT_NAME}
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/test_solvers
    COMMAND ${CMAKE_COMMAND} -E echo "Cleaned build directory"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Removing all build artifacts"
)