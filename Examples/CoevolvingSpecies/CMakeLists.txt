# Set the minimum version of CMake
cmake_minimum_required(VERSION 3.13)

# Set the project's name
project(CoevolvingSpecies LANGUAGES CXX CUDA)

# Set cuda, curand paths
set(CURAND_DIR /usr/local/cuda/lib64)
set(CUDA_DIR /usr/local/cuda/lib64)

# Import code base
if(NOT SPEGG_ROOT)
	set(SPEGG_ROOT /home/wormlab/OneDrive/Personal_Project/spegg	)
	list(
		APPEND 
		CMAKE_MODULE_PATH 
		${SPEGG_ROOT}
	)
endif()

# Add cmake
include(${SPEGG_ROOT}/FindSPEGG.cmake)

# Find CUDA and SPEGG packages
find_package(CUDA REQUIRED)
find_package(SPEGG REQUIRED)

# Set cuda architecture and optimization flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 --extended-lambda --expt-relaxed-constexpr")
set(CMAKE_CUDA_ARCHITECTURES 50 61 70 75 80 86)

# Set where the binary would go
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

# Add include dirs
include_directories(${CMAKE_SOURCE_DIR}/header)
include_directories(${SPEGG_INCLUDE_DIRS})
include_directories(${CUDA_INCLUDE_DIRS})

# Define src directory
set(SRC ${CMAKE_SOURCE_DIR}/src)

# Add executables
add_executable(
	${PROJECT_NAME} 
	${SRC}/coevolutionSimulator.cu 
	${SRC}/coevolvingSpeciesParents.cu 
	${SRC}/myConsumer_genotype_phenotype_maps.cu 
	${SRC}/myConsumer.cu 
	${SRC}/deme_specific_data_class_coevolvingSpecies.cu 
	${SRC}/coevolvingSpecie.cu 
	${SRC}/myResource.cu 
	${SRC}/UpdateBehavior.cu 
	${SRC}/myResource_genotype_phenotype_maps.cu 
	${SRC}/migrate_coevolvingSpecies.cu 
	${SRC}/genotype_phenotype_map.cu 
	${SRC}/update_coevolvingSpecies.cu 
	${SRC}/MigrationBehavior.cu 
	main.cpp
)

# Enable separable compilation for all CUDA targets
set_target_properties(
	${PROJECT_NAME} 
	PROPERTIES 
	CUDA_SEPARABLE_COMPILATION 
	ON
)

# Target link to the executable
target_link_libraries(
	${PROJECT_NAME} 
	PRIVATE 
	SPEGG::spegg_codebase
)

