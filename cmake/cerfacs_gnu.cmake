####################################################################
# COMPILER
####################################################################

set ( CMAKE_C_COMPILER       mpicc  )
set ( CMAKE_CXX_COMPILER     mpicxx )
set ( CMAKE_Fortran_COMPILER mpif90 )

####################################################################
# OpenMP FLAGS
####################################################################

set( OpenMP_C_FLAGS             "-fopenmp" )
set( OpenMP_CXX_FLAGS           "-fopenmp" )
set( OpenMP_Fortran_FLAGS       "-fopenmp" )

####################################################################
# FYPP PREPROCESSOR FOR INTEL (line numbering warning)
####################################################################

set( FYPP_NO_LINE_NUMBERING TRUE )

####################################################################
# COMMON FLAGS
####################################################################

set( ECBUILD_C_FLAGS        "${CMAKE_C_FLAGS} -Wall -Wno-deprecated-declarations -Wno-parentheses" )
set( ECBUILD_CXX_FLAGS      "${CMAKE_CXX_FLAGS} -Wall -Wno-deprecated-declarations -Wno-parentheses" )
set( ECBUILD_Fortran_FLAGS  "")

####################################################################
# RELEASE FLAGS
####################################################################

set( ECBUILD_C_FLAGS_RELEASE        "-O3 -g -DNDEBUG" )
set( ECBUILD_CXX_FLAGS_RELEASE      "-O3 -g -DNDEBUG" )
set( ECBUILD_Fortran_FLAGS_RELEASE  "-O3 -g -funroll-all-loops -finline-functions -DNDEBUG" )

####################################################################
# RELEASE WITH DEBUG INFO FLAGS
####################################################################

set( ECBUILD_C_FLAGS_RELWITHDEBINFO        "-O3 -g -fbacktrace" )
set( ECBUILD_CXX_FLAGS_RELWITHDEBINFO      "-O3 -g -fbacktrace" )
set( ECBUILD_Fortran_FLAGS_RELWITHDEBINFO  "-O3 -g -fbacktrace -funroll-all-loops -finline-functions" )

####################################################################
# BIT (REPRODUCIBLE) FLAGS
####################################################################

set( ECBUILD_C_FLAGS_BIT       "-O2" )
set( ECBUILD_CXX_FLAGS_BIT     "-O2" )
set( ECBUILD_Fortran_FLAGS_BIT "-O2 -funroll-all-loops -finline-functions" )

####################################################################
# PRODUCTION FLAGS
####################################################################

set( ECBUILD_C_FLAGS_PRODUCTION       "-O3 -DNDEBUG" )
set( ECBUILD_CXX_FLAGS_PRODUCTION     "-O3 -DNDEBUG" )
set( ECBUILD_Fortran_FLAGS_PRODUCTION "-O3 -funroll-all-loops -finline-functions -DNDEBUG" )

####################################################################
# DEBUG FLAGS
####################################################################

set( ECBUILD_C_FLAGS_DEBUG        "-O0 -g -fbacktrace" )
set( ECBUILD_CXX_FLAGS_DEBUG      "-O0 -g -fbacktrace" )
set( ECBUILD_Fortran_FLAGS_DEBUG  "-O0 -g -Wextra -Wall -ftrapv -fall-intrinsics -fcheck=all -fimplicit-none -ffpe-trap=invalid,zero,overflow,denormal -fcheck-array-temporaries -finit-derived -finit-integer=-999 -finit-real=snan -fbacktrace" )

####################################################################
# AVOID THE DETECTION OF SYSTEM BOOST
####################################################################

set( Boost_NO_BOOST_CMAKE on )
set( Boost_NO_SYSTEM_PATHS on )

####################################################################
# LINK FLAGS
####################################################################

set( ECBUILD_Fortran_LINK_FLAGS  "-lblas" )
set( ECBUILD_CXX_LINK_FLAGS "-lblas")

###################################################################
#
# PACKAGES
###################################################################

set( NETCDF_PATH  $ENV{NETCDF4_DIR})
set( HDF5_PATH    $ENV{HDF5_DIR})
