/*
 * (C) Crown copyright 2021, Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_OBSAUXPARAMETERS_H_
#define AQ_OBSAUXPARAMETERS_H_

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

namespace aq {

/// Parameters taken by the ObsAuxControl, ObsAuxControlCorrection and ObsAuxCovariance classes.

// -----------------------------------------------------------------------------

class ObsAuxCovarianceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsAuxCovarianceParameters, Parameters)

 public:
  oops::OptionalParameter<double> insitu{"insitu", this};
};

// -----------------------------------------------------------------------------

class ObsAuxControlParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsAuxControlParameters, Parameters)

 public:
  oops::OptionalParameter<double> insitu{"insitu", this};
  oops::OptionalParameter<ObsAuxCovarianceParameters> covariance{"covariance", this};
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_OBSAUXPARAMETERS_H_
