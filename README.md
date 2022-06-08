# AQ
Air Quality 4DEnVar suite for the SEEDS project based on the JEDI Data Assimilation environment

# Obtaining AQ
<pre>
git clone git@github.com:andreapiacentini/aq.git
</pre>

# Installation
In order to compile AQ and its prerequisites from JSCDA

## Set the following environment variables
`AQ_JEDI_SRC` - source directory of aq (the one where the git clone has been extracted)
`AQ_JEDI_BLD` - ecbuild directory of the aq bundle
`AQ_JEDI_DIR` - installation directory of the aq and prerequisite include, libs, exe
`AQ_JEDI_UTI` - installation directory of compilation utils
`AQ_JEDI_TEST TIER` - (optional) level of testing, default to 1. 2 for 3d tests
`AQ_JEDI_LOCAL_PATH_TESTFILES` - (optional) path of test data on the local machine

## Load the appropriate modules
On kraken:
<pre>
Currently Loaded Modulefiles:
  1) compiler/gcc/8.3.0
  2) compiler/intel/19.0.5.281
  3) mpi/intelmpi/2019.4.243
  4) lib/phdf5/1.10.4_impi
  5) lib/netcdf-fortran/4.4.4_phdf5_1.10.4
  6) tools/cmake/3.17.0(default)
  7) lib/boost/1.66.0_intel
  8) python/3.7.7
</pre>

On nemo with intel compilers:
<pre>
Currently Loaded Modulefiles:
  1) compiler/intel/18.0.1.163
  2) mpi/intelmpi/2018.1.163-ddt
  3) lib/phdf5/1.10.4_impi
  4) lib/netcdf-fortran/4.4.4_phdf5_1.10.4
  5) compiler/gcc/5.4.0
  6) tools/cmake/3.23.3
  7) lib/boost/1.66.0_gcc540
  8) python/2.7-shared
  9) python/3.7.7
</pre>

On nemo with gnu compilers:
<pre>
Currently Loaded Modulefiles:
  1) compiler/gcc/9.4.0
  2) mpi/intelmpi/2021.5.1
  3) lib/phdf5/1.10.8_gcc940_impi
  4) lib/netcdf-fortran/4.4.4_phdf5_1.10.8_gcc940
  5) tools/cmake/3.23.2
  6) lib/boost/1.66.0_gcc940
  7) python/3.7.7
</pre>

On davinci:
<pre>
Currently Loaded Modulefiles:
  1) compiler/intel/19.0.5.281
  2) mpi/intelmpi/2018.1.163
  3) lib/phdf5/1.10.4_impi
  4) lib/netcdf-fortran/4.4.4_phdf5_1.10.4
  5) lib/boost/1.66.0_intel
  6) python/3.9.0
</pre>

## Recover ecbuild
<pre>
mkdir -p $AQ_JEDI_UTI
cd $AQ_JEDI_UTI
git clone https://github.com/ecmwf/ecbuild.git
</pre>

## Set PATH for using ecbuild
<pre>
export PATH=${AQ_JEDI_UTI}/ecbuild/bin:$PATH
</pre>

You can already set the path for the utilities that will be compiled together with the bundle
<pre>
export PATH=${AQ_JEDI_DIR}/bin:$PATH
</pre>

## Define extra paths in the .bashrc
`NetCDF_INCLUDE_DIRS` is undefined and required by `ecbuild:ecbuild_find_package`
<pre>
export NetCDF_INCLUDE_DIRS=$NETCDF_DIR/include
</pre>or
<pre>
export NetCDF_INCLUDE_DIRS=$NETCDF4_INCDIR
</pre>

Make BOOST point to the appropriate recent version (`module load lib/boost/1.66.0_intel`)
<pre>
export BOOST_ROOT=/softs/local_intel/boost/1.66.0
</pre>or on nemo
<pre>
export BOOST_ROOT=/data/softs/local/boost/1.66.0_gcc540
</pre>

On nemo make Eigen3 point to the appropriate recent version
<pre>
export Eigen3_ROOT=/data/softs/local/eigen/3.3.5_gcc540
</pre>

## Point to the appropriate toolchain for compilation options
For the intel based compilations
<pre>
export ECBUILD_TOOLCHAIN_DIR=$AQ_JEDI_SRC/cmake
export ECBUILD_TOOLCHAIN=${ECBUILD_TOOLCHAIN_DIR}/cerfacs_intel.cmake
</pre>

For the gnu based compilations
<pre>
export ECBUILD_TOOLCHAIN_DIR=$AQ_JEDI_SRC/cmake
export ECBUILD_TOOLCHAIN=${ECBUILD_TOOLCHAIN_DIR}/cerfacs_gnu.cmake
</pre>

## Install the bundle by executing the install script
From `$AQ_JEDI_SRC`
<pre>
./install_aq.ksh
</pre>

## Run the tests
From `$AQ_JEDI_BLD`
<pre>
cd aq
ctest
</pre>
