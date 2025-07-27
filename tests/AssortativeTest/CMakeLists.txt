cmake_minimum_required(VERSION 3.13)
project(AssortativeTest LANGUAGES CXX CUDA)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

# Set SPEGG location (like GTEST_ROOT)
if (NOT SPEGG_ROOT)
    set(SPEGG_ROOT "~/Documents/DatNewProject/spegg")
    # Add custom module path
    list(APPEND CMAKE_MODULE_PATH "${SPEGG_ROOT}")
endif()

if (NOT GTEST_ROOT)
    set(GTEST_ROOT "~/googletest")
endif()

find_package(Threads REQUIRED)
find_package(GTest REQUIRED)
find_package(SPEGG REQUIRED)
# Find CUDA package
find_package(CUDA REQUIRED)

# Set Required C++ Standard Ensure the appropriate C++ standard is set explicitly
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set C and C++ compilers explicitly
set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

# Set CUDA compiler explicitly
set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")

# Set cuda architecture and optimization flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 --extended-lambda --expt-relaxed-constexpr")
set(CMAKE_CUDA_ARCHITECTURES 61-real)

include_directories(${GTEST_INCLUDE_DIRS})
include_directories(${GTEST_ROOT}/googlemock/include)
include_directories(${SPEGG_INCLUDE_DIRS})  # SPEGG headers
# Add include directories for headers
include_directories(${CMAKE_SOURCE_DIR})

# Add CUDA include directories explicitly
include_directories(/usr/local/cuda/include)
include_directories(${CUDA_INCLUDE_DIRS})
# Resolve external libraries
# Find the libraries
find_library(LIBCONFIG_LIB config++ REQUIRED)
find_library(CURAND_LIB curand REQUIRED HINTS /usr/local/cuda/lib64)
find_library(CUDA_LIB cuda REQUIRED HINTS /usr/local/cuda/lib64)
find_library(RT_LIB rt REQUIRED)

add_executable(${PROJECT_NAME} ${CODEBASE_LIBS} GenPhenMapTest.h GenPhenMapTest.cu AssortativeNeonatesTest.cu)

# Enable separable compilation for all CUDA targets
set_target_properties(${PROJECT_NAME} PROPERTIES 
    CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(${PROJECT_NAME} PRIVATE GTest::gmock GTest::GTest Threads::Threads SPEGG::spegg_codebase  
    ${CODEBASE_LIBS}
    ${LIBCONFIG_LIB}
    ${CURAND_LIB}
    ${CUDA_LIB}
    ${RT_LIB})
