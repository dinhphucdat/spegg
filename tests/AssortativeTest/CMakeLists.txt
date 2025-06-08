cmake_minimum_required(VERSION 3.18)  # Higher version for better CUDA support
project(AssortativeTest LANGUAGES CXX CUDA)

# Include required modules
include(FetchContent)

# Set output directory (fixed syntax)
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/../")

# Set Required C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set CUDA standard
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Set compilers explicitly (if needed)
# Note: Usually better to let CMake find them automatically
# set(CMAKE_C_COMPILER "/usr/bin/gcc")
# set(CMAKE_CXX_COMPILER "/usr/bin/g++")

# CUDA architecture and optimization flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 --extended-lambda --expt-relaxed-constexpr")
set(CMAKE_CUDA_ARCHITECTURES "61-real")

# Find packages
find_package(Threads REQUIRED)
find_package(CUDA REQUIRED)

# Fetch GoogleTest if not provided
if (NOT GTEST_ROOT)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0  # Use a more recent stable version
    )
    FetchContent_MakeAvailable(googletest)
    set(GTEST_FOUND ON)
else()
    find_package(GTest REQUIRED)
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/../../header)
include_directories(${CMAKE_SOURCE_DIR})

# Add CUDA include directories
include_directories(${CUDA_INCLUDE_DIRS})

# Find external libraries
find_library(LIBCONFIG_LIB config++ REQUIRED)
find_library(CURAND_LIB curand REQUIRED HINTS ${CUDA_TOOLKIT_ROOT_DIR}/lib64)
find_library(CUDART_LIB cudart REQUIRED HINTS ${CUDA_TOOLKIT_ROOT_DIR}/lib64)
find_library(RT_LIB rt REQUIRED)

# Find codebase libraries
set(CODEBASE_DIR ${CMAKE_SOURCE_DIR}/../../objdir2)
file(GLOB CODEBASE_LIBS "${CODEBASE_DIR}/*.a")

# Verify that codebase libraries were found
if(NOT CODEBASE_LIBS)
    message(WARNING "No codebase libraries found in ${CODEBASE_DIR}")
endif()

# Add library search path (should come before add_executable)
link_directories(${CODEBASE_DIR})

# Define source files
set(SOURCES
    GenPhenMapTest.cu
    AssortativeNeonatesTest.cu
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Enable separable compilation for CUDA
set_target_properties(${PROJECT_NAME} PROPERTIES 
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE
    ${CODEBASE_LIBS}
    ${LIBCONFIG_LIB}
    ${CURAND_LIB}
    ${CUDART_LIB}
    ${RT_LIB}
    Threads::Threads
)

# Link GoogleTest libraries
if(GTEST_FOUND)
    if(TARGET GTest::gtest)
        target_link_libraries(${PROJECT_NAME} PRIVATE GTest::gtest GTest::gmock)
    else()
        target_link_libraries(${PROJECT_NAME} PRIVATE gtest gmock)
    endif()
endif()

# Compiler options (uncomment if needed)
# target_compile_options(${PROJECT_NAME} PRIVATE 
#     $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra>
#     $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-Wall,-Wextra>
# )

# Debug information
message(STATUS "CUDA Version: ${CUDA_VERSION}")
message(STATUS "CUDA Architectures: ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "Found ${CMAKE_CUDA_COMPILER_ID} CUDA compiler: ${CMAKE_CUDA_COMPILER}")
message(STATUS "Codebase libraries: ${CODEBASE_LIBS}")
