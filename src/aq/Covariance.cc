/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "aq/Covariance.h"

#include "aq/Geometry.h"
#include "aq/Increment.h"
#include "aq/State.h"
#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace aq {
// -----------------------------------------------------------------------------
Covariance::Covariance(const Geometry & resol, const oops::Variables & vars,
                                     const eckit::Configuration & conf,
                                     const State &, const State &) {
  oops::Log::error() << "Covariance: constructor not implemented" << std::endl;
}
// -----------------------------------------------------------------------------
Covariance::~Covariance() {
  oops::Log::error() << "Covariance: destructor not implemented" << std::endl;
}
// -----------------------------------------------------------------------------
void Covariance::multiply(const Increment & dxin, Increment & dxout) const {
  oops::Log::error() << "Covariance: multiply not implemented" << std::endl;
}
// -----------------------------------------------------------------------------
void Covariance::inverseMultiply(const Increment & dxin, Increment & dxout) const {
  oops::Log::error() << "Covariance: inverseMultiply not implemented" << std::endl;
}
// -----------------------------------------------------------------------------
void Covariance::randomize(Increment & dx) const {
  oops::Log::error() << "Covariance: randomize not implemented" << std::endl;
}
// -----------------------------------------------------------------------------
void Covariance::print(std::ostream & os) const {
  os << "Covariance::print not implemented";
}
// -----------------------------------------------------------------------------

}  // namespace aq
