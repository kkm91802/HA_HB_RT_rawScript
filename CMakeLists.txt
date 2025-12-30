cmake_minimum_required(VERSION 3.10)

project(shared_folder C)
set(CMAKE_C_STANDARD 99)

# Collect all source files in src/
file(GLOB SOURCES "src/*.c")

# Build the executable
add_executable(shared_folder ${SOURCES})

# Include headers
target_include_directories(shared_folder PRIVATE ${PROJECT_SOURCE_DIR}/include)

# ------------------------------------------------------
# Find and link FFTW3 + math library (consistent style)
# ------------------------------------------------------
find_library(FFTW_LIB fftw3 REQUIRED)
if(FFTW_LIB)
    message(STATUS "Found FFTW3 library: ${FFTW_LIB}")
    target_link_libraries(shared_folder PRIVATE ${FFTW_LIB})
else()
    message(FATAL_ERROR "FFTW3 library not found! Install with: sudo apt install libfftw3-dev")
endif()

# Always link math library too
target_link_libraries(shared_folder PRIVATE m)

# ------------------------------------------------------
# Optional: Python interpreter
# ------------------------------------------------------
find_package(Python3 COMPONENTS Interpreter REQUIRED)
