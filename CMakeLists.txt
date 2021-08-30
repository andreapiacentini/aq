# (C) Copyright 2021 Cerfacs and UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# ####################
# aq data assimilation
# ####################

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )
project( aq VERSION 1.0.0 LANGUAGES CXX Fortran )

## Ecbuild integration
find_package( ecbuild 3.5 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild)
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

# Default release mode
set( ECBUILD_DEFAULT_BUILD_TYPE Release )

## Configuration options
option( ENABLE_AQ_DOC "Build AQ documentation" OFF )
option( ENABLE_MKL "Use MKL for LAPACK implementation (if available)" ON )
option( ENABLE_GPTL "Use GPTL profiling library (if available)" OFF )
option( ENABLE_AUTOPROFILING "Enable function-based autoprofiling with GPTL (if available)" OFF )

include( ${PROJECT_NAME}_compiler_flags )
include( GNUInstallDirs )

## Dependencies
find_package( jedicmake QUIET )  # Prefer find modules from jedi-cmake

if( ENABLE_MKL )
    find_package( MKL )
endif()
if( MKL_FOUND )
    set( LAPACK_LIBRARIES ${MKL_LIBRARIES} )
else()
    find_package( LAPACK REQUIRED )
endif()

find_package( Eigen3 REQUIRED NO_MODULE HINTS
    $ENV{Eigen3_ROOT} $ENV{EIGEN3_ROOT} $ENV{Eigen_ROOT} $ENV{EIGEN_ROOT}
    $ENV{Eigen3_PATH} $ENV{EIGEN3_PATH} $ENV{Eigen_PATH} $ENV{EIGEN_PATH})

find_package( OpenMP COMPONENTS CXX Fortran )
find_package( MPI REQUIRED COMPONENTS C CXX Fortran )
find_package( NetCDF REQUIRED COMPONENTS Fortran )
find_package( HDF5 REQUIRED COMPONENTS HL Fortran )
find_package( Boost 1.64.0 REQUIRED)
find_package( eckit 1.11.6 REQUIRED COMPONENTS MPI )
find_package( fckit 0.7.0 REQUIRED )
if(OpenMP_FOUND)
    find_package( atlas 0.20.2 REQUIRED COMPONENTS OMP OMP_Fortran )
else()
    find_package( atlas 0.20.2 REQUIRED )
endif()

#Optional Dependencies
if( jedicmake_FOUND AND ENABLE_GPTL )
    find_package( GPTL )
    add_compile_definitions(ENABLE_GPTL)
endif()

if ( GPTL_FOUND AND ENABLE_AUTOPROFILING )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_AUTOPROFILING_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Fortran_AUTOPROFILING_FLAGS}")
endif()

# JSON validation
find_package(nlohmann_json QUIET)
find_package(nlohmann_json_schema_validator QUIET)

## Sources
add_subdirectory( src )
add_subdirectory( test )


## Package Config

ecbuild_install_project( NAME ${PROJECT_NAME} )

#Export Targets
set(EXPORTED_TARGET_NAMESPACE "") #TODO: set to "aq::" once all packages are ready for namespaced targets
set(INSTALL_TARGETS_PATH lib/cmake/${PROJECT_NAME})
get_filename_component(TARGETS_FILE "${PROJECT_TARGETS_FILE}" NAME)

#Export package config for build and install trees
set( BINDIR_BUILD_EXPORT ${CMAKE_BINARY_DIR}/bin ) #Exported binary (tools) directory locations for build-tree
file( MAKE_DIRECTORY ${BINDIR_BUILD_EXPORT} )
set( BINDIR_INSTALL_EXPORT ${CMAKE_INSTALL_FULL_BINDIR} ) #Exported binary (tools) directory locations for install-tree
ecbuild_generate_project_config(${PROJECT_NAME}-config.cmake.in PATH_VARS BINDIR_BUILD_EXPORT BINDIR_INSTALL_EXPORT )

ecbuild_print_summary()
