# (C) Copyright 2017-2020 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# #########
# aq bundle
# #########

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

find_package( ecbuild 3.5 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild)

project( aq-bundle VERSION 1.0.0 LANGUAGES C CXX Fortran )

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include( ecbuild_bundle )

# Default release mode
set( ECBUILD_DEFAULT_BUILD_TYPE Release )

# Enable MPI
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI")

# Initialize
ecbuild_bundle_initialize()
ecbuild_bundle( PROJECT jedicmake GIT "https://github.com/jcsda/jedi-cmake.git" BRANCH develop UPDATE )
include( jedicmake/cmake/Functions/git_functions.cmake  )

# ECMWF libs
# ----------
option("BUNDLE_SKIP_ECKIT" "Don't build eckit" "OFF" ) # Skip eckit build unless user passes -DBUNDLE_SKIP_ECKIT=OFF
option("BUNDLE_SKIP_FCKIT" "Don't build fckit" "OFF") # Skip fckit build unless user passes -DBUNDLE_SKIP_FCKIT=OFF
option("BUNDLE_SKIP_ATLAS" "Don't build atlas" "OFF") # Skip atlas build unless user passes -DBUNDLE_SKIP_ATLAS=OFF

ecbuild_bundle( PROJECT eckit GIT "https://github.com/ecmwf/eckit.git" TAG 1.16.0 )
ecbuild_bundle( PROJECT fckit GIT "https://github.com/ecmwf/fckit.git" TAG 0.9.2 )
ecbuild_bundle( PROJECT atlas GIT "https://github.com/ecmwf/atlas.git" TAG 0.24.1 )

# Core JEDI repositories
# ----------------------
ecbuild_bundle( PROJECT oops  GIT "https://github.com/jcsda/oops.git"  BRANCH develop UPDATE )
#ecbuild_bundle( PROJECT oops  GIT "https://github.com/andreapiacentini/oops.git"  BRANCH bugfix/atlas_test UPDATE )
ecbuild_bundle( PROJECT saber GIT "https://github.com/jcsda/saber.git" BRANCH develop UPDATE )
#ecbuild_bundle( PROJECT ioda  GIT "https://github.com/jcsda/ioda.git"  BRANCH develop UPDATE )
#ecbuild_bundle( PROJECT ufo   GIT "https://github.com/jcsda/ufo.git"   BRANCH develop UPDATE )

# AQ repository
# ------------------------------------
ecbuild_bundle( PROJECT aq SOURCE .. )

# ioda, ufo and saber test data
#----------------------------------------

# If IODA branch is being built set GIT_BRANCH_FUNC to IODA's current branch.
# If a tagged version of IODA is being built set GIT_TAG_FUNC to ioda's current tag. In this case,
# IODA test files will be download from UCAR DASH and ioda-data repo will not be cloned.
# When LOCAL_PATH_JEDI_TESTFILES is set to the directory of IODA test files stored
# in a local directory, ioda-data repo will not be cloned

#find_branch_name(REPO_DIR_NAME ioda)
# When LOCAL_PATH_JEDI_TESTFILES is set to the directory of IODA test files stored
# in a local directory, ioda-data repo will not be cloned
#if( NOT DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES} AND NOT DEFINED ${GIT_TAG_FUNC} )
#  ecbuild_bundle( PROJECT ioda-data GIT "https://github.com/JCSDA-internal/ioda-data.git" BRANCH develop UPDATE )
#endif()

# If IODA's current branch is available in ioda-data repo, that branch will be checked out
#branch_checkout (REPO_DIR_NAME ioda-data
#                 BRANCH ${GIT_BRANCH_FUNC} )

# same procedure for ufo-data
#find_branch_name(REPO_DIR_NAME ufo)
#if( NOT DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES} AND NOT DEFINED ${GIT_TAG_FUNC} )
#  ecbuild_bundle( PROJECT ufo-data GIT "https://github.com/JCSDA-internal/ufo-data.git" BRANCH develop UPDATE )
#endif()

# If UFO's current branch is available in ufo-data repo, that branch will be checked out
#branch_checkout (REPO_DIR_NAME ufo-data
#                 BRANCH ${GIT_BRANCH_FUNC} )

# same procedure for saber-data
find_branch_name(REPO_DIR_NAME saber)
if( NOT DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES} AND NOT DEFINED ${GIT_TAG_FUNC} )
  ecbuild_bundle( PROJECT saber-data GIT "https://github.com/JCSDA-internal/saber-data.git" BRANCH develop UPDATE )
endif()

# If saber's current branch is available in saber-data repo, that branch will be checked out
branch_checkout (REPO_DIR_NAME saber-data
                 BRANCH ${GIT_BRANCH_FUNC} )

# Build Doxygen documentation
# ---------------------------
option(BUILD_AQ_BUNDLE_DOC "Build documentation" OFF)
if(BUILD_AQ_BUNDLE_DOC)
   add_subdirectory( Documentation )
endif(BUILD_AQ_BUNDLE_DOC)

ecbuild_bundle_finalize()