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
