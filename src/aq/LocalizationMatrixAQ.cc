/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/LocalizationMatrixAQ.h"

#include "eckit/config/Configuration.h"
#include "model/GeometryAQ.h"
#include "model/IncrementAQ.h"
// AQ #include "model/AqFortran.h"

// -----------------------------------------------------------------------------
namespace aq {
// -----------------------------------------------------------------------------
LocalizationMatrixAQ::LocalizationMatrixAQ(const GeometryAQ & resol,
                                           const eckit::Configuration & config) {
}
// -----------------------------------------------------------------------------
LocalizationMatrixAQ::~LocalizationMatrixAQ() {
}
// -----------------------------------------------------------------------------
void LocalizationMatrixAQ::randomize(IncrementAQ & dx) const {
}
// -----------------------------------------------------------------------------
void LocalizationMatrixAQ::multiply(IncrementAQ & dx) const {
  IncrementAQ dxtmp(dx);
}
// -----------------------------------------------------------------------------
void LocalizationMatrixAQ::print(std::ostream & os) const {
  os << "LocalizationMatrixAQ::print not implemented";
}
// -----------------------------------------------------------------------------

}  // namespace aq
