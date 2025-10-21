# Either spegg location or spegg root is okay, but the latter is preferred
if (SPEGG_LOCATION AND NOT SPEGG_ROOT)
    set(SPEGG_ROOT ${SPEGG_LOCATION})
endif()

# Expand tilde in root
if (SPEGG_ROOT)
    get_filename_component(SPEGG_ROOT ${SPEGG_ROOT} ABSOLUTE)
endif()

# Find the root directory
find_path(SPEGG_ROOT_DIR 
    NAMES CMakeLists.txt 
    HINTS ${SPEGG_ROOT} $ENV{SPEGG_ROOT} 
    PATH_SUFFIXES . 
    DOC "Spegg root directory."
)

add_subdirectory(
    ../../pybind11
    ${CMAKE_BINARY_DIR}/pybind11${CMAKE_BINARY_DIR}/pybind11-build
)

if (SPEGG_ROOT_DIR)
    # Find include directories
    set(SPEGG_INCLUDE_DIRS
        ${SPEGG_ROOT_DIR}/header
        ${SPEGG_ROOT_DIR}/header/environ
        ${SPEGG_ROOT_DIR}/header/math
        ${SPEGG_ROOT_DIR}/header/species
        ${SPEGG_ROOT_DIR}/header/species/add_kids
        ${SPEGG_ROOT_DIR}/header/species/movement
        ${SPEGG_ROOT_DIR}/header/species/update
        ${SPEGG_ROOT_DIR}/header/util
    )
    # Find the compiled libraries directory
    set(SPEGG_LIB_DIR ${SPEGG_ROOT_DIR}/lib)
    set(CURAND_DIR "/usr/local/cuda/lib64")
    # Find all static libraries
    file(GLOB SPEGG_LIBRARIES "${SPEGG_LIB_DIR}/*.a")
    # Check if all binary files are found
    if (SPEGG_LIBRARIES)
        set(SPEGG_FOUND TRUE)
        # Find external dependencies that SPEGG needs
        find_library(SPEGG_LIBCONFIG_LIB config++ REQUIRED)
        message(STATUS "curand dir is ${CURAND_DIR}")
        find_library(SPEGG_CURAND_LIB curand REQUIRED HINTS ${CURAND_DIR})
        find_library(SPEGG_CUDA_LIB cuda REQUIRED HINTS ${CUDA_DIR})
        find_library(SPEGG_RT_LIB rt REQUIRED)
        list(
            APPEND PYBIND11_LIB 
            pybind11::headers
            pybind11::embed
            pybind11::module)
        
        # Add external libraries to SPEGG_LIBRARIES
        list(APPEND SPEGG_LIBRARIES
            ${SPEGG_LIBCONFIG_LIB}
            ${SPEGG_CURAND_LIB}
            ${SPEGG_CUDA_LIB}
            ${SPEGG_RT_LIB}
            ${PYBIND11_LIB}
        )
        
        # Create imported target if not already created
        if(NOT TARGET SPEGG::spegg_codebase)
            add_library(SPEGG::spegg_codebase INTERFACE IMPORTED)
            set_target_properties(SPEGG::spegg_codebase PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${SPEGG_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${SPEGG_LIBRARIES}"
            )
        endif()
    else()
        set(SPEGG_FOUND FALSE)
    endif()
else()
    set(SPEGG_FOUND FALSE)
endif()

# Handle find_package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SPEGG
    FOUND_VAR SPEGG_FOUND
    REQUIRED_VARS SPEGG_ROOT_DIR SPEGG_LIBRARIES SPEGG_INCLUDE_DIRS
    FAIL_MESSAGE "Could not find SPEGG library. Set SPEGG_ROOT or SPEGG_LOCATION to the SPEGG installation directory."
)

# Mark variables as advanced
mark_as_advanced(
    SPEGG_ROOT_DIR
    SPEGG_INCLUDE_DIRS
    SPEGG_LIBRARIES
    SPEGG_LIBCONFIG_LIB
    SPEGG_CURAND_LIB
    SPEGG_CUDA_LIB
    SPEGG_RT_LIB
    PYBIND11_LIB
)
