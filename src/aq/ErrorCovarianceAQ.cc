/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "aq/ErrorCovarianceAQ.h"

#include "aq/GeometryAQ.h"
#include "aq/IncrementAQ.h"
#include "aq/StateAQ.h"
#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace aq {
// -----------------------------------------------------------------------------
ErrorCovarianceAQ::ErrorCovarianceAQ(const GeometryAQ & resol, const oops::Variables & vars,
                                     const eckit::Configuration & conf,
                                     const StateAQ &, const StateAQ &) {
  oops::Log::error() << "ErrorCovarianceAQ: constructor not implemented" << std::endl;
}
// -----------------------------------------------------------------------------
ErrorCovarianceAQ::~ErrorCovarianceAQ() {
  oops::Log::error() << "ErrorCovarianceAQ: destructor not implemented" << std::endl;
}
// -----------------------------------------------------------------------------
void ErrorCovarianceAQ::multiply(const IncrementAQ & dxin, IncrementAQ & dxout) const {
  oops::Log::error() << "ErrorCovarianceAQ: multiply not implemented" << std::endl;
}
// -----------------------------------------------------------------------------
void ErrorCovarianceAQ::inverseMultiply(const IncrementAQ & dxin, IncrementAQ & dxout) const {
  oops::Log::error() << "ErrorCovarianceAQ: inverseMultiply not implemented" << std::endl;
}
// -----------------------------------------------------------------------------
void ErrorCovarianceAQ::randomize(IncrementAQ & dx) const {
  oops::Log::error() << "ErrorCovarianceAQ: randomize not implemented" << std::endl;
}
// -----------------------------------------------------------------------------
void ErrorCovarianceAQ::print(std::ostream & os) const {
  os << "ErrorCovarianceAQ::print not implemented";
}
// -----------------------------------------------------------------------------

}  // namespace aq
