/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
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
