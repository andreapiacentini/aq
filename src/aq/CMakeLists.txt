include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
list( APPEND aq_src_files
  Traits.h
  ChangeVar.cc
  ChangeVar.h
  LinearChangeVar.cc
  LinearChangeVar.h
  Covariance.cc
  Covariance.h
  Fields.cc
  Fields.h
  Geometry.cc
  Geometry.h
  GeometryIterator.cc
  GeometryIterator.h
  GetValues.cc
  GetValues.h
  LinearGetValues.cc
  LinearGetValues.h
  GeoVals.cc
  GeoVals.h
  Increment.cc
  Increment.h
  Locations.cc
  Locations.h
  ModelAuxControl.h
  ModelAuxCovariance.h
  ModelAuxIncrement.h
  #AQ Do we really need them ?
  ObsAuxControl.cc
  ObsAuxControl.h
  ObsAuxParameters.h
  ObsAuxCovariance.cc
  ObsAuxCovariance.h
  ObsAuxIncrement.cc
  ObsAuxIncrement.h
  ObsAuxPreconditioner.cc
  ObsAuxPreconditioner.h
  #AQ end of question
  ObsData.h
  ObsDiagnostics.h
  ObsIterator.cc
  ObsIterator.h
  ObsOpBase.cc
  ObsOpBase.h
  LinearObsOpBase.cc
  LinearObsOpBase.h
  ObsOperator.cc
  ObsOperator.h
  LinearObsOperator.cc
  LinearObsOperator.h
  ObsSpace.cc
  ObsSpace.h
  ObsInSitu.cc
  ObsInSitu.h
  LinearObsInSitu.cc
  LinearObsInSitu.h
  ObsVec.cc
  ObsVec.h
  State.cc
  State.h
  instantiateAqChangeVarFactory.h
  aq_blas_mod.F90
  aq_constants_mod.F90
  aq_fields_interface.F90
  aq_fields_mod.F90
  aq_field_io_mod.F90
  aq_geom_interface.F90
  aq_geom_mod.F90
  aq_geom_iter_interface.F90
  aq_geom_iter_mod.F90
  aq_getvalues_interface.F90
  aq_getvalues_mod.F90
  aq_geovals_interface.F90
  aq_geovals_mod.F90
  aq_interp_mod.F90
  aq_locs_f.cc
  aq_locs_f.h
  aq_locs_mod.F90
  aq_obsdb_interface.F90
  aq_obsdb_mod.F90
  aq_obsvec_interface.F90
  aq_obsvec_mod.F90
  aq_obsvec_random_f.h
  aq_obsvec_random_f.cc
  aq_insitu_interface.F90
  aq_insitu_mod.F90
  aq_tools_mod.F90
  AnalyticInit.cc
  AnalyticInit.h
  FinalCheck.cc
  FinalCheck.h
  h5_read_mod.F90
  h5_selection_mod.F90
  h5_utils_mod.F90
  h5_write_mod.F90
)

# Requirement sources
set( AQ_SRC_DEP ${aq_src_files} )

# Requirement libraries
set( AQ_LIB_DEP oops saber atlas atlas_f ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} )

# Optional OpenMP
if( OpenMP_Fortran_FOUND )
  list( APPEND AQ_LIB_DEP OpenMP::OpenMP_Fortran )
endif()

# Requirement includes
set( AQ_INC_DEP ${HDF5_INCLUDE_DIRS} )

# Add the library
ecbuild_add_library( TARGET          aq
                     SOURCES         ${AQ_SRC_DEP}
                     PUBLIC_LIBS     ${AQ_LIB_DEP}
                     PUBLIC_INCLUDES ${AQ_INC_DEP}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${AQ_LINKER_LANGUAGE}
                   )

target_compile_features( aq PUBLIC cxx_std_11 )
target_include_directories( aq PUBLIC "$<BUILD_INTERFACE:${AQ_INCLUDE_DIRS}>" "$<BUILD_INTERFACE:${AQ_EXTRA_INCLUDE_DIRS}>" )
