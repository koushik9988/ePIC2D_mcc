cmake_minimum_required(VERSION 3.22)
project(ePIC++ LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -O2")

# Source files for main executable
file(GLOB SOURCES "src/*.cpp" "cross_section/*.cpp")

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

# Find OpenMP (optional but recommended for performance)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "Found OpenMP: ${OpenMP_CXX_VERSION}")
else()
    message(WARNING "OpenMP not found. Building without OpenMP support.")
endif()

# Find FFTW3
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED IMPORTED_TARGET fftw3)

if(NOT FFTW_FOUND)
    message(FATAL_ERROR "FFTW3 not found. Please install FFTW3 (e.g., sudo apt install libfftw3-dev)")
endif()

message(STATUS "Found FFTW3: ${FFTW_VERSION}")

# Print version information
message(STATUS "Found Python version: ${Python3_VERSION}")
message(STATUS "Python executable: ${Python3_EXECUTABLE}")
message(STATUS "Python include dirs: ${Python3_INCLUDE_DIRS}")
message(STATUS "NumPy include dirs: ${Python3_NumPy_INCLUDE_DIRS}")
message(STATUS "Found HDF5 version: ${HDF5_VERSION}")
message(STATUS "HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")


# Include directories
include_directories(
    include
    linearalgebra
    cross_section
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
    PkgConfig::FFTW
)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()
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
    PkgConfig::FFTW
)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_solvers PRIVATE OpenMP::OpenMP_CXX)
endif()
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
