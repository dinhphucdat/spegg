# Set the minimum version of CMake
cmake_minimum_required(VERSION 3.10)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Force fPIC for standard C++ files
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
# Set cuda architecture and optimization flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -O3 -fPIC --extended-lambda --expt-relaxed-constexpr")
set(CMAKE_CUDA_ARCHITECTURES 50 61 70 75 80 86)

# Set the project name
project(spegg LANGUAGES CXX CUDA)

########### CONSIDER EDITTING THESE TWO VARIABLES! #################
set(CURAND_DIR_CACHE /usr/local/cuda/lib64 CACHE PATH "Path to the installed curand library")
set(CUDA_DIR_CACHE /usr/local/cuda/lib64 CACHE PATH "Path to the installed cuda library")
###################################################################
set(CURAND_DIR ${CURAND_DIR_CACHE})
set(CUDA_DIR ${CUDA_DIR_CACHE})

# Find CUDA package
find_package(CUDA REQUIRED)

# Find pybind11
if (NOT TARGET pybind11::module)
    add_subdirectory(
        ${CMAKE_SOURCE_DIR}/pybind11 
        ${CMAKE_CURRENT_BINARY_DIR}/pybind11_base
    )
endif()

# Set Required C++ Standard Ensure the appropriate C++ standard is set explicitly
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set C and C++ compilers explicitly
set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

# Set CUDA compiler explicitly
set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")

# Set where all the output would go
set(FULL_OUTPUT "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${FULL_OUTPUT}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${FULL_OUTPUT}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${FULL_OUTPUT}/..")

# Add include directories for headers
include_directories(${CMAKE_SOURCE_DIR}/header)
include_directories(${CMAKE_SOURCE_DIR}/header/environ)
include_directories(${CMAKE_SOURCE_DIR}/header/math)
include_directories(${CMAKE_SOURCE_DIR}/header/species)
include_directories(${CMAKE_SOURCE_DIR}/header/species/add_kids)
include_directories(${CMAKE_SOURCE_DIR}/header/species/movement)
include_directories(${CMAKE_SOURCE_DIR}/header/species/update)
include_directories(${CMAKE_SOURCE_DIR}/header/util)

# Add CUDA include directories explicitly
include_directories(/usr/local/cuda/include)
include_directories(${CUDA_INCLUDE_DIRS})

# Define source directories
set(SRC ${CMAKE_SOURCE_DIR}/src)
set(CUDA ${SRC}/cuda)
set(ENVIRON ${CUDA}/environ)
set(MATH ${CUDA}/math)
set(SPECIES ${CUDA}/species)
set(UTIL ${CUDA}/util)
set(ADD_KIDS ${SPECIES}/add_kids)
set(MOVEMENT ${SPECIES}/movement)
set(UPDATE ${SPECIES}/update)

add_library(spegg_obj STATIC 
    ${SRC}/Simulation_Class.cu ${SPECIES}/inds.cu ${SPECIES}/inds_stochastic.cu 
    ${SPECIES}/inds_stochastic_migratory.cu ${ENVIRON}/environment.cu 
    ${MATH}/mating_thrust_prob_table_demes.cu ${MATH}/thrust_prob_table_demes.cu 
    ${MATH}/thrust_prob_table.cu ${MATH}/histogram.cu ${UTIL}/footimer2.cpp 
    ${UTIL}/amplify.cu ${UTIL}/reduce_by_key_with_zeroes.cu ${UTIL}/Shuffle_by_Key.cu 
    ${MATH}/random_variables_functions.cu ${MATH}/statistics_class.cu 
    ${MATH}/demographic_statistics_class.cu ${SPECIES}/deme_specific_data_class.cu 
    ${ADD_KIDS}/parents_class.cu ${ADD_KIDS}/neonates_class.cu 
    ${ADD_KIDS}/genotype_phenotype_map_parameters.cu ${ADD_KIDS}/genetic_deme_specific_data.cu 
    ${ADD_KIDS}/assortative_mating_neonates_class.cu ${ADD_KIDS}/assortative_mating_parents_class.cu 
    ${UPDATE}/determine_mortality.cu ${MOVEMENT}/MigrationFunctions.cu 
    ${UTIL}/one_dim_two_dim.cu ${UTIL}/remove_duplicate_pairs.cu ${UTIL}/which_function.cu 
    ${UTIL}/gather_values_by_deme.cu ${UTIL}/Sampling_Input.cu ${UTIL}/Sampling_Event.cu 
    ${UTIL}/Sample_without_Replacement_1Pass.cu ${UTIL}/Sample_without_Replacement_byDeme.cu 
    ${UTIL}/Sample_With_Replacement.cu ${UTIL}/Sample_without_Replacement_Test.cu 
    ${ADD_KIDS}/species_specific_mate_sampling_rules.cu
)

# This is the "Double-Lock"
target_compile_options(spegg_obj PRIVATE 
    $<$<COMPILE_LANGUAGE:CXX>:-fPIC>
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler -fPIC>
)

# Enable separable compilation for all CUDA targets
set_target_properties(spegg_obj PROPERTIES 
    CUDA_SEPARABLE_COMPILATION ON 
    CUDA_RESOLVE_DEVICE_SYMBOLS ON 
    POSITION_INDEPENDENT_CODE ON 
    CXX_VISIBILITY_PRESET default
    CUDA_VISIBILITY_PRESET default
    VISIBILITY_INLINES_HIDDEN OFF
    CUDA_RUNTIME_LIBRARY Shared
)

# Resolve external libraries
# Find the libraries
find_library(LIBCONFIG_LIB config++ REQUIRED)
find_library(CURAND_LIB curand REQUIRED HINTS ${CURAND_DIR})
find_library(CUDA_LIB cuda REQUIRED HINTS ${CUDA_DIR})
find_library(RT_LIB rt REQUIRED)

# Pybind11 libraries
list(
    APPEND PYBIND11_LIB 
    pybind11::headers
    pybind11::embed
    pybind11::module)

# Store the list of external libraries in a variable
set(EXTERNAL_LIBS
    ${LIBCONFIG_LIB}
    ${CURAND_LIB}
    ${CUDA_LIB}
    ${RT_LIB}
    ${PYBIND11_LIB}
    -Wl,--whole-archive 
    -Wl,--no-whole-archive
    cudart 
    cudadevrt
)

# Link target libraries
target_link_libraries(spegg_obj ${EXTERNAL_LIBS})
# Set the appropriate compile options for CUDA files
# set_target_properties(a.out PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
# Add a custom clean target to remove the object directory and executable
# By renaming the target to `clean_custom`, you avoid conflicting with the reserved `clean` target.
add_custom_target(clean_obj
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${FULL_OUTPUT}
    COMMENT "Cleaning object directory"
)
#     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}