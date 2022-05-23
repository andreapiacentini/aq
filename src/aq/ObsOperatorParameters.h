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

#ifndef AQ_OBSOPERATORPARAMETERS_H_
#define AQ_OBSOPERATORPARAMETERS_H_

#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/Parameters.h"

namespace aq {

/// Parameters controlling the observation operator for the AQ model.
class ObservationParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObservationParameters, Parameters);
 public:
  oops::ConfigurationParameter config{this};
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_OBSOPERATORPARAMETERS_H_
