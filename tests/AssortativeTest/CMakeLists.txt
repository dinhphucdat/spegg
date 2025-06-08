cmake_minimum_required(VERSION 3.18)
project(AssortativeTest LANGUAGES CXX CUDA)

# Include required modules
include(FetchContent)

# Set standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# CUDA settings
set(CMAKE_CUDA_ARCHITECTURES "61")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 --extended-lambda --expt-relaxed-constexpr")

# Find packages
find_package(Threads REQUIRED)

# Fetch GoogleTest
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/../../header)
include_directories(${CMAKE_SOURCE_DIR})

# Find libraries (only essential ones for now)
find_library(LIBCONFIG_LIB config++)
find_library(CURAND_LIB curand HINTS /usr/local/cuda/lib64)

# Find codebase libraries
set(CODEBASE_DIR ${CMAKE_SOURCE_DIR}/../../objdir2)
file(GLOB CODEBASE_LIBS "${CODEBASE_DIR}/*.a")

# Define sources
set(SOURCES
    GenPhenMapTest.cu
    AssortativeNeonatesTest.cu
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Set CUDA properties (minimal approach)
set_target_properties(${PROJECT_NAME} PROPERTIES 
    CUDA_SEPARABLE_COMPILATION ON  # Try without separable compilation first
)

# Link only essential libraries
target_link_libraries(${PROJECT_NAME} PRIVATE
    gtest
    gmock
    ${CODEBASE_LIBS}
)

# Add other libraries if the basic build works
if(LIBCONFIG_LIB)
    target_link_libraries(${PROJECT_NAME} PRIVATE ${LIBCONFIG_LIB})
endif()

if(CURAND_LIB)
    target_link_libraries(${PROJECT_NAME} PRIVATE ${CURAND_LIB})
endif()

# Debug output
message(STATUS "CUDA Architectures: ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "Codebase libraries found: ${CODEBASE_LIBS}")