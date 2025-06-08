cmake_minimum_required(VERSION 3.13)
project(AssortativeTest LANGUAGES CXX CUDA)

set(${EXECUTABLE_OUTPUT_PATH} "${CMAKE_SOURCE_DIR}../")

# if (NOT GTEST_ROOT)
#     set(GTEST_ROOT "~/googletest/googletest")
# endif()


# Fetch GoogleTest if not provided
if (NOT GTEST_ROOT)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1 # Replace with your desired version
    )
    FetchContent_MakeAvailable(googletest)
endif()

find_package(Threads REQUIRED)
find_package(GTest REQUIRED)
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
# set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")

# Set cuda architecture and optimization flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 --extended-lambda --expt-relaxed-constexpr")
set(CMAKE_CUDA_ARCHITECTURES 61-real)

include_directories(${GTEST_INCLUDE_DIRS})

# Add include directories for headers
include_directories(${CMAKE_SOURCE_DIR}/../../header)
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

# Find shared object files from the codebase
set(CODEBASE ${CMAKE_SOURCE_DIR}/../../objdir2)
file(GLOB CODEBASE_LIBS "${CODEBASE}/*.a")

add_executable(${PROJECT_NAME} ${CODEBASE_LIBS} GenPhenMapTest.h GenPhenMapTest.cu AssortativeNeonatesTest.cu)

# Enable separable compilation for all CUDA targets
set_target_properties(${PROJECT_NAME} PROPERTIES 
    CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(${PROJECT_NAME} PRIVATE GTest::GTest Threads::Threads ${CODEBASE_LIBS} 
    ${CODEBASE_LIBS}
    ${LIBCONFIG_LIB}
    ${CURAND_LIB}
    ${CUDA_LIB}
    ${RT_LIB})

# Add library search path
link_directories(${CODEBASE})

# target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Werror)
