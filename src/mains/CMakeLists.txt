ecbuild_add_executable( TARGET  aq_hofx3d.x
                        SOURCES aqHofX3D.cc
                        LIBS    aq
                      )

ecbuild_add_executable( TARGET  aq_diffstates.x
                        SOURCES aqDiffStates.cc
                        LIBS    aq
                      )

ecbuild_add_executable( TARGET  aq_addincrement.x
                        SOURCES aqAddIncrement.cc
                        LIBS    aq
                      )

ecbuild_add_executable( TARGET  aq_ens_mean_variance.x
                        SOURCES aqEnsMeanVariance.cc
                        LIBS    aq
                      )

ecbuild_add_executable( TARGET  aq_ens_recenter.x
                        SOURCES aqEnsRecenter.cc
                        LIBS    aq
                      )

ecbuild_add_executable( TARGET  aq_error_covariance_toolbox.x
                        SOURCES aqErrorCovarianceToolbox.cc
                        LIBS    aq
                                saber
                      )

ecbuild_add_executable( TARGET  aq_var.x
                        SOURCES aqVar.cc
                        LIBS    aq
                                saber
                      )
