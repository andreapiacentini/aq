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

#include <cmath>

#include "aq/FieldsAQ.h"
#include "aq/GeometryAQ.h"
#include "aq/IncrementAQ.h"
#include "aq/StateAQ.h"
#include "eckit/config/Configuration.h"
#include "oops/assimilation/GMRESR.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace aq {
// -----------------------------------------------------------------------------
ErrorCovarianceAQ::ErrorCovarianceAQ(const GeometryAQ & resol, const oops::Variables & vars,
                                     const eckit::Configuration & conf,
                                     const StateAQ &, const StateAQ &) {
  oops::Log::trace() << "ErrorCovarianceAQ created" << std::endl;
}
// -----------------------------------------------------------------------------
ErrorCovarianceAQ::~ErrorCovarianceAQ() {
  oops::Log::trace() << "ErrorCovarianceAQ destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ErrorCovarianceAQ::multiply(const IncrementAQ & dxin, IncrementAQ & dxout) const {
}
// -----------------------------------------------------------------------------
void ErrorCovarianceAQ::inverseMultiply(const IncrementAQ & dxin, IncrementAQ & dxout) const {
  oops::IdentityMatrix<IncrementAQ> Id;
  dxout.zero();
  GMRESR(dxout, dxin, *this, Id, 20, 1.0e-5);
}
// -----------------------------------------------------------------------------
void ErrorCovarianceAQ::randomize(IncrementAQ & dx) const {
}
// -----------------------------------------------------------------------------
void ErrorCovarianceAQ::print(std::ostream & os) const {
  os << "ErrorCovarianceAQ::print not implemented";
}
// -----------------------------------------------------------------------------

}  // namespace aq
