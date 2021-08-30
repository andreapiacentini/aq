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
LocalizationMatrixAQ.cc
LocalizationMatrixAQ.h
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
instantiateAqLocalizationFactory.h
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

ecbuild_add_library( TARGET   aq
                     SOURCES  ${aq_src_files}
                     PUBLIC_LIBS     oops
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/aq
                     INSTALL_HEADERS LISTED
                     PUBLIC_INCLUDES ${HDF5_INCLUDE_DIRS}
                     LINKER_LANGUAGE CXX )

#Configure include directory layout for build-tree to match install-tree
set(AQ_BUILD_DIR_INCLUDE_PATH ${CMAKE_BINARY_DIR}/${PROJECT_NAME}/include)
add_custom_target(aq_headers ALL COMMAND ${CMAKE_COMMAND} -E make_directory "${AQ_BUILD_DIR_INCLUDE_PATH}/oops"
                                    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR} ${AQ_BUILD_DIR_INCLUDE_PATH}/oops/aq)
add_dependencies(aq aq_headers)
target_include_directories(aq PUBLIC $<BUILD_INTERFACE:${AQ_BUILD_DIR_INCLUDE_PATH}>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../..>)

#Fortran module output directory for build interface
set(AQ_MODULE_DIR module/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION}/aq)
set_target_properties(aq PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${AQ_MODULE_DIR})
install(DIRECTORY ${CMAKE_BINARY_DIR}/${AQ_MODULE_DIR}/ DESTINATION ${AQ_MODULE_DIR})
target_include_directories(aq INTERFACE
                                $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${AQ_MODULE_DIR}>
                                $<INSTALL_INTERFACE:${AQ_MODULE_DIR}>)
