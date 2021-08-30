list( APPEND aq_data_files
  E0001/HMGEMS02+2010070401.nc
  E0002/HMGEMS02+2010070401.nc
  E0003/HMGEMS02+2010070401.nc
  E0004/HMGEMS02+2010070401.nc
  E0005/HMGEMS02+2010070401.nc
  E0006/HMGEMS02+2010070401.nc
  E0007/HMGEMS02+2010070401.nc
  E0008/HMGEMS02+2010070401.nc
  E0009/HMGEMS02+2010070401.nc
  E0010/HMGEMS02+2010070401.nc
  HDAT+20100704.h5
  HGOM+20100704_CO.gom.nc
  HGOM+20100704_O3.gom.nc
)

list( APPEND aq_testinput
  #AQ testinput/3densvar.yaml
  #AQ testinput/4densvar.yaml
  testinput/addincrement.yaml
  testinput/addincrement_scaled.yaml
  testinput/diffstates.yaml
  testinput/dirac_loc_3d.yaml
  #AQ testinput/dirac_loc_4d.yaml
  testinput/dirac_no_loc.yaml
  testinput/ens_recenter.yaml
  testinput/ens_variance.yaml
  testinput/ens_variance_inflation_value.yaml
  testinput/getvalues.yaml
  testinput/hofx3d.yaml
  testinput/interfaces.yaml
  testinput/lineargetvalues.yaml
  testinput/parameters_loc_3d.yaml
)

list( APPEND aq_testoutput
  #AQ testoutput/3densvar.test
  #AQ testoutput/4densvar.test
  testoutput/addincrement.test
  testoutput/addincrement_scaled.test
  testoutput/diffstates.test
  testoutput/dirac_loc_3d.test
  #AQ testoutput/dirac_loc_4d.test
  testoutput/dirac_no_loc.test
  testoutput/ens_recenter.test
  testoutput/ens_variance.test
  testoutput/ens_variance_inflation_value.test
  testoutput/hofx3d.test
  testoutput/parameters_loc_3d.test
)

# Link scripts
list( APPEND tools_files
    aq_data_downloader.sh
)

foreach( FILENAME ${tools_files} )
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                             ${PROJECT_SOURCE_DIR}/test/${FILENAME}
                             ${CMAKE_BINARY_DIR}/bin/${FILENAME} )
endforeach()


# Create data directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

# Create data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${aq_testinput})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create data directory for reference output and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
foreach(FILENAME ${aq_testoutput})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

ecbuild_add_resources( TARGET   aq_test_scripts
                       SOURCES_PACK
                       ${aq_testinput}
                     )

#####################################################################
# download test
#####################################################################

ecbuild_add_test( TARGET    get_data_aq
                  TYPE      SCRIPT
                  COMMAND   ${CMAKE_BINARY_DIR}/bin/aq_data_downloader.sh
                  ARGS  ${aq_data_files} )

#####################################################################
# interface tests
#####################################################################

ecbuild_add_test( TARGET  test_aq_geometry
                  MPI 6
                  OMP 4
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    aq )

ecbuild_add_test( TARGET  test_aq_geovals
                  MPI 6
                  OMP 4
                  SOURCES executables/TestGeoVaLs.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    aq )

ecbuild_add_test( TARGET  test_aq_state
                  MPI 6
                  OMP 4
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    aq )

ecbuild_add_test( TARGET  test_aq_increment
                  MPI 6
                  OMP 4
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    aq )

ecbuild_add_test( TARGET test_aq_locations
                  MPI 6
                  OMP 4
                  SOURCES executables/TestLocations.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    aq )

ecbuild_add_test( TARGET test_aq_obsspace
                  MPI 6
                  OMP 4
                  SOURCES executables/TestObsSpace.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    aq )

ecbuild_add_test( TARGET test_aq_obs_vector
                  MPI 6
                  OMP 4
                  SOURCES executables/TestObsVector.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    aq )

ecbuild_add_test( TARGET test_aq_obsdatavector
                  MPI 6
                  OMP 4
                  SOURCES executables/TestObsDataVector.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    aq )

ecbuild_add_test( TARGET  test_aq_obs_operator
                  MPI 6
                  OMP 4
                  SOURCES executables/TestObsOperator.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    aq )

ecbuild_add_test( TARGET  test_aq_lin_obs_operator
                  MPI 6
                  OMP 4
                  SOURCES executables/TestLinearObsOperator.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    aq )

ecbuild_add_test( TARGET test_aq_obs_error_covar
                  MPI 6
                  OMP 4
                  SOURCES executables/TestObsErrorCovariance.cc
                  ARGS "testinput/interfaces.yaml"
                  LIBS    aq )

ecbuild_add_test( TARGET  test_aq_obs_aux
                  MPI 6
                  OMP 4
                  SOURCES executables/TestObsAuxControl.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    aq )

ecbuild_add_test( TARGET  test_aq_obs_aux_increment
                  MPI 6
                  OMP 4
                  SOURCES executables/TestObsAuxIncrement.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    aq )

ecbuild_add_test( TARGET  test_aq_obs_aux_covariance
                  MPI 6
                  OMP 4
                  SOURCES executables/TestObsAuxCovariance.cc
                  ARGS    "testinput/interfaces.yaml"
                  LIBS    aq )

ecbuild_add_test( TARGET  test_aq_getvalues
                  MPI 6
                  OMP 4
                  SOURCES executables/TestGetValues.cc
                  ARGS    "testinput/getvalues.yaml"
                  LIBS    aq )

ecbuild_add_test( TARGET  test_aq_lineargetvalues
                  MPI 6
                  OMP 4
                  SOURCES executables/TestLinearGetValues.cc
                  ARGS    "testinput/lineargetvalues.yaml"
                  LIBS    aq )

#####################################################################
# obs-related tests
#####################################################################

ecbuild_add_test( TARGET test_aq_hofx3d
                  MPI 6
                  OMP 4
                  ARGS testinput/hofx3d.yaml
                  COMMAND  aq_hofx3d.x )

#####################################################################
# ensemble-related tests
#####################################################################

ecbuild_add_test( TARGET test_aq_ens_variance
                  MPI 6
                  OMP 4
                  ARGS testinput/ens_variance.yaml
                  COMMAND  aq_ens_variance.x )

ecbuild_add_test( TARGET test_aq_ens_recenter
                  MPI 6
                  OMP 4
                  ARGS testinput/ens_recenter.yaml
                  COMMAND  aq_ens_recenter.x )

ecbuild_add_test( TARGET test_aq_ens_variance_inflation_value
                  MPI 6
                  OMP 4
                  ARGS testinput/ens_variance_inflation_value.yaml
                  COMMAND  aq_ens_variance.x )

#####################################################################
# B Matrix tests
#####################################################################

ecbuild_add_test( TARGET test_aq_parameters_loc_3d
                  MPI 6
                  OMP 4
                  ARGS testinput/parameters_loc_3d.yaml
                  COMMAND  aq_parameters.x )

#####################################################################
# dirac tests
#####################################################################

ecbuild_add_test( TARGET test_aq_dirac_loc_3d
                  MPI 6
                  OMP 4
                  ARGS testinput/dirac_loc_3d.yaml
                  COMMAND  aq_dirac.x )

ecbuild_add_test( TARGET test_aq_dirac_no_loc
                  MPI 6
                  OMP 4
                  ARGS testinput/dirac_no_loc.yaml
                  COMMAND  aq_dirac.x)

#####################################################################
# 3d variational tests
#####################################################################

#ecbuild_add_test( TARGET test_aq_3densvar
#                  OMP 2
#                  ARGS testinput/3densvar.yaml
#                  COMMAND  aq_4dvar.x )

######################################################################
## 4d variational tests
######################################################################

#ecbuild_add_test( TARGET test_aq_4densvar
#                  MPI 7
#                  ARGS testinput/4densvar.yaml
#                  COMMAND  aq_4dvar.x )

######################################################################
## state-related tests
######################################################################

ecbuild_add_test( TARGET test_aq_diffstates
                  MPI 6
                  OMP 4
                  ARGS testinput/diffstates.yaml
                  COMMAND  aq_diffstates.x
		  )

ecbuild_add_test( TARGET test_aq_addincrement
                  MPI 6
                  OMP 4
                  ARGS testinput/addincrement.yaml
                  COMMAND  aq_addincrement.x
                  TEST_DEPENDS test_aq_diffstates )

ecbuild_add_test( TARGET test_aq_addincrement_scaled
                  MPI 6
                  OMP 4
                  ARGS testinput/addincrement_scaled.yaml
                  COMMAND  aq_addincrement.x
                  TEST_DEPENDS test_aq_diffstates )
