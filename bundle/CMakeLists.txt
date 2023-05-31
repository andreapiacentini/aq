# (C) Copyright 2017-2020 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# #########
# aq bundle
# #########

cmake_minimum_required( VERSION 3.12 )
project( aq-bundle VERSION 0.0.1 LANGUAGES C CXX Fortran )

find_package(ecbuild)
set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_CURRENT_SOURCE_DIR}/ecbuild/cmake;${CMAKE_MODULE_PATH}")
include( ecbuild_bundle )

# Default release mode
set( ECBUILD_DEFAULT_BUILD_TYPE RelWithDebInfo )

# Enable MPI
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

# Avoid rebuild after each commit
set( ECBUILD_RECORD_GIT_COMMIT_SHA1 OFF CACHE BOOL
     "When building ecbuild projects that are Git repos, create variables recording the full and short Git revision" )

# Initialize
ecbuild_bundle_initialize()
ecbuild_bundle( PROJECT jedicmake GIT "https://github.com/jcsda/jedi-cmake.git" BRANCH develop UPDATE RECURSIVE )
include( jedicmake/cmake/Functions/git_functions.cmake  )

# ECMWF libs
# ----------
option("BUNDLE_SKIP_ECKIT" "Don't build eckit" "OFF" ) # Skip eckit build unless user passes -DBUNDLE_SKIP_ECKIT=OFF
option("BUNDLE_SKIP_FCKIT" "Don't build fckit" "OFF") # Skip fckit build unless user passes -DBUNDLE_SKIP_FCKIT=OFF
option("BUNDLE_SKIP_ATLAS" "Don't build atlas" "OFF") # Skip atlas build unless user passes -DBUNDLE_SKIP_ATLAS=OFF

ecbuild_bundle( PROJECT eckit GIT "https://github.com/ecmwf/eckit.git" TAG 1.20.2 )
ecbuild_bundle( PROJECT fckit GIT "https://github.com/ecmwf/fckit.git" TAG 0.9.5 )
ecbuild_bundle( PROJECT atlas GIT "https://github.com/ecmwf/atlas.git" TAG 0.31.1 )

# Core JEDI repositories
# ----------------------
ecbuild_bundle( PROJECT oops  GIT "git@github.com:JCSDA-internal/oops.git"  BRANCH feature/Bapps UPDATE )
ecbuild_bundle( PROJECT saber GIT "git@github.com:JCSDA-internal/saber.git" BRANCH feature/Bapps UPDATE )

# AQ repository
# ------------------------------------
ecbuild_bundle( PROJECT aq SOURCE .. )

ecbuild_bundle_finalize()
