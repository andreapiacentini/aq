include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
list( APPEND aq_src_files
  AqTraits.h
  # The ChangeVar (direct and TLAD are not needed for AQ
  # but they could be required by generic oops tests)
  #AQ ChangeVarAQ.cc
  #AQ ChangeVarAQ.h
  #AQ ChangeVarTLADAQ.cc
  #AQ ChangeVarTLADAQ.h
  ErrorCovarianceAQ.cc
  ErrorCovarianceAQ.h
  FieldsAQ.cc
  FieldsAQ.h
  GeometryAQ.cc
  GeometryAQ.h
  GeometryAQIterator.cc
  GeometryAQIterator.h
  GetValuesAQ.cc
  GetValuesAQ.h
  GetValuesTLAD.cc
  GetValuesTLAD.h
  GomAQ.cc
  GomAQ.h
  IncrementAQ.cc
  IncrementAQ.h
  LocationsAQ.cc
  LocationsAQ.h
  #AQ Do we really need them ?
  ObsBias.cc
  ObsBias.h
  ObsBiasParameters.h
  ObsBiasCovariance.cc
  ObsBiasCovariance.h
  ObsBiasIncrement.cc
  ObsBiasIncrement.h
  #AQ end of question
  ObsDataAQ.h
  ObsDiagsAQ.h
  ObsIteratorAQ.cc
  ObsIteratorAQ.h
  ObsOpBaseAQ.cc
  ObsOpBaseAQ.h
  ObsOpBaseTLAD.cc
  ObsOpBaseTLAD.h
  ObsOperatorAQ.cc
  ObsOperatorAQ.h
  ObsOperatorTLAD.cc
  ObsOperatorTLAD.h
  ObsSpaceAQ.cc
  ObsSpaceAQ.h
  ObsStreamAQ.cc
  ObsStreamAQ.h
  ObsStreamTLAD.cc
  ObsStreamTLAD.h
  ObsVecAQ.cc
  ObsVecAQ.h
  StateAQ.cc
  StateAQ.h
  #AQ instantiateAqChangeVarFactory.h
  AqFortran.h
  aq_blas_mod.F90
  #AQ aq_change_var_interface.F90
  #AQ aq_change_var_mod.F90
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
  aq_gom_interface.F90
  aq_gom_mod.F90
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
  aq_stream_interface.F90
  aq_stream_mod.F90
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


# Requirement libraries
set(AQ_LIB_DEP oops saber atlas atlas_f ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} )

# Optional OpenMP
if(OpenMP_Fortran_FOUND)
  list( APPEND AQ_LIB_DEP OpenMP::OpenMP_Fortran)
endif()

# Requirement sources
set(AQ_SRC_DEP ${aq_src_files} )


# Add the library
ecbuild_add_library( TARGET          aq
                     SOURCES         ${AQ_SRC_DEP}
                     PUBLIC_LIBS     ${AQ_LIB_DEP}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${AQ_LINKER_LANGUAGE}
                   )

target_compile_features( aq PUBLIC cxx_std_11 )
target_include_directories( aq PUBLIC "$<BUILD_INTERFACE:${AQ_INCLUDE_DIRS}>" "$<BUILD_INTERFACE:${AQ_EXTRA_INCLUDE_DIRS}>" )
