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

ecbuild_add_executable( TARGET  aq_dirac.x
                        SOURCES aqDirac.cc
                        LIBS    aq
                                saber
                      )

ecbuild_add_executable( TARGET  aq_ens_variance.x
                        SOURCES aqEnsVariance.cc
                        LIBS    aq
                      )

ecbuild_add_executable( TARGET  aq_ens_recenter.x
                        SOURCES aqEnsRecenter.cc
                        LIBS    aq
                      )

ecbuild_add_executable( TARGET  aq_error_covariance_training.x
                        SOURCES aqErrorCovarianceTraining.cc
                        LIBS    aq
                                saber
                      )

ecbuild_add_executable( TARGET  aq_var.x
                        SOURCES aqVar.cc
                        LIBS    aq
                                saber
                      )
