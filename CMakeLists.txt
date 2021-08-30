# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# ####################
# aq data assimilation
# ####################

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

project( aq VERSION 0.0.1 LANGUAGES C CXX Fortran )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH})

set( CMAKE_DIRECTORY_LABELS "aq" )

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_OS_TESTS           OFF CACHE BOOL "Disable OS tests" FORCE )
set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.5 )

################################################################################
# Project
################################################################################

ecbuild_declare_project()

ecbuild_enable_fortran( REQUIRED )

set( AQ_LINKER_LANGUAGE CXX )

################################################################################
# Standard compiler flags
################################################################################

include( aq_extra_macros )
option( ENABLE_AQ_DOC "Build AQ documentation" OFF )

################################################################################
# Dependencies
################################################################################

# jedi-cmake
find_package( jedicmake QUIET )

# HDF5
find_package( HDF5 REQUIRED COMPONENTS Fortran_HL HL Fortran )

# NetCDF
find_package( NetCDF REQUIRED COMPONENTS Fortran )

# atlas
find_package( atlas 0.20.2 REQUIRED )

# oops
find_package( oops 1.0.0 REQUIRED )

# saber
find_package( saber 1.0.0 REQUIRED )

# OpenMP
find_package( OpenMP COMPONENTS Fortran )

# MPI
ecbuild_find_mpi( COMPONENTS CXX Fortran REQUIRED )

################################################################################
# Export package info
################################################################################

list( APPEND AQ_TPLS atlas atlas_f oops ufo saber )

set( AQ_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src
                          ${CMAKE_CURRENT_BINARY_DIR}/src
                          ${CMAKE_Fortran_MODULE_DIRECTORY})

set( AQ_LIBRARIES aq )

get_directory_property( AQ_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${AQ_TPLS} )
  string( TOUPPER ${_tpl} TPL )
  list( APPEND AQ_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
  list( APPEND AQ_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
  list( APPEND AQ_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

################################################################################
# Sources
################################################################################

add_subdirectory( src )
add_subdirectory( test )

if( ENABLE_AQ_DOC )
    add_subdirectory( Documents )
endif()

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME aq )

# print the summary of the configuration
ecbuild_print_summary()
