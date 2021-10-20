####################################################################
# ARCHITECTURE
####################################################################

set( Fortran_AUTOMATIC_ARRAYS_LIMIT 32768 )  # (32 kb)
math( EXPR Fortran_AUTOMATIC_ARRAYS_LIMIT_KB "${Fortran_AUTOMATIC_ARRAYS_LIMIT}/1024" )
set( Fortran_FLAG_STACK_ARRAYS     "-no-heap-arrays" )
set( Fortran_FLAG_AUTOMATIC_ARRAYS "-heap-arrays ${Fortran_AUTOMATIC_ARRAYS_LIMIT_KB}" )

####################################################################
# COMPILER
####################################################################

set ( CMAKE_C_COMPILER       mpiicc  )
set ( CMAKE_CXX_COMPILER     mpiicpc )
set ( CMAKE_Fortran_COMPILER mpiifort)

####################################################################
# OpenMP FLAGS
####################################################################

set( OpenMP_C_FLAGS             "-qopenmp" )
set( OpenMP_CXX_FLAGS           "-qopenmp" )
set( OpenMP_Fortran_FLAGS       "-qopenmp" )

####################################################################
# FYPP PREPROCESSOR FOR INTEL (line numbering warning)
####################################################################

set( FYPP_NO_LINE_NUMBERING TRUE )

####################################################################
# COMMON FLAGS
####################################################################

# for diagnostics:
#  -diag-enable=vec -diag-file -Winline

set( ECBUILD_C_FLAGS        "-xAVX" )
set( ECBUILD_CXX_FLAGS      "-xAVX -std=c++14" )
set( ECBUILD_Fortran_FLAGS  "-xAVX -r8 -mkl ${Fortran_FLAG_AUTOMATIC_ARRAYS}")

####################################################################
# RELEASE FLAGS
####################################################################

# for diagnostics:
#  -diag-enable=vec -diag-file -Winline

set( ECBUILD_C_FLAGS_RELEASE        "-O3 -g -DNDEBUG -finline-limit=500" )
set( ECBUILD_CXX_FLAGS_RELEASE      "-O3 -g -DNDEBUG -finline-limit=500 -std=c++14" )
set( ECBUILD_Fortran_FLAGS_RELEASE  "-O3 -g -DNDEBUG -nowarn -unroll -inline -finline-limit=500 -align array64byte" )

####################################################################
# RELEASE WITH DEBUG INFO FLAGS
####################################################################

# for diagnostics:
#  -diag-enable=vec -diag-file -Winline

set( ECBUILD_C_FLAGS_RELWITHDEBINFO        "-O3 -g -traceback -finline-limit=500" )
set( ECBUILD_CXX_FLAGS_RELWITHDEBINFO      "-O3 -g -traceback -finline-limit=500 -std=c++14" )
set( ECBUILD_Fortran_FLAGS_RELWITHDEBINFO  "-O3 -g -traceback -nowarn -unroll -inline -finline-limit=500 -align array64byte" )

####################################################################
# BIT (REPRODUCIBLE) FLAGS
####################################################################

set( ECBUILD_C_FLAGS_BIT       "-O2 -g -DNDEBUG -fp-speculation=strict -fp-model precise" )
set( ECBUILD_CXX_FLAGS_BIT     "-O2 -g -DNDEBUG -fp-speculation=strict -fp-model precise -std=c++14" )
set( ECBUILD_Fortran_FLAGS_BIT "-O2 -g -DNDEBUG -nowarn -unroll -inline -fp-speculation=strict -fp-model precise" )

####################################################################
# PRODUCTION FLAGS
####################################################################

set( ECBUILD_C_FLAGS_PRODUCTION       "-O3 -g -fp-model precise" )
set( ECBUILD_CXX_FLAGS_PRODUCTION     "-O3 -g -fp-model precise -std=c++14" )
set( ECBUILD_Fortran_FLAGS_PRODUCTION "-O3 -g -nowarn -fp-model precise" )

####################################################################
# DEBUG FLAGS
####################################################################

set( ECBUILD_C_FLAGS_DEBUG        "-O0 -g -traceback -fp-trap=common" )
set( ECBUILD_CXX_FLAGS_DEBUG      "-O0 -g -traceback -fp-trap=common -std=c++14" )
# -check all implies -check bounds
set( ECBUILD_Fortran_FLAGS_DEBUG  "-O0 -g -traceback -fp-model precise -warn unused -warn all -heap-arrays -fpe0 -ftrapuv -check all" )
#set( ECBUILD_Fortran_FLAGS_DEBUG  "-O0 -g -traceback -fp-model precise -heap-arrays -fpe0 -ftrapuv -check noarg_temp_created" )

####################################################################
# AVOID THE DETECTION OF SYSTEM BOOST
####################################################################

set( Boost_NO_BOOST_CMAKE on )
set( Boost_NO_SYSTEM_PATHS on )

####################################################################
# LINK FLAGS
####################################################################

set( ECBUILD_Fortran_LINK_FLAGS  "-mkl" )

###################################################################
# 
# PACKAGES
###################################################################

set( NETCDF_PATH  $ENV{NETCDF4_DIR})
set( HDF5_PATH    $ENV{HDF5_DIR})
set( LAPACK_FOUND $ENV{MKLROOT})
