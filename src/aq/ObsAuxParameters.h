/*
 * This file is part of the Air Quality Ensemble Data Assimilation suite AQ.
 *
 * (C) Copyright 2022 CERFACS
 *
 * AQ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * AQ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * A copy of the GNU Lesser General Public License is distributed
 * along with AQ (files LICENSE.md, COPYING and COPYING.LESSER).
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
