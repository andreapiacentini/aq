/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "aq/ObsOperatorAQ.h"

#include "aq/GomAQ.h"
#include "aq/LocationsAQ.h"
#include "aq/ObsBias.h"
#include "aq/ObsDiagsAQ.h"
#include "aq/ObsOpBaseAQ.h"
#include "aq/ObsSpaceAQ.h"
#include "aq/ObsVecAQ.h"
#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"

namespace aq {

// -----------------------------------------------------------------------------

ObsOperatorAQ::ObsOperatorAQ(const ObsSpaceAQ & os, const eckit::Configuration & conf)
  : oper_(ObsOpFactory::create(os, conf))
{}

// -----------------------------------------------------------------------------

ObsOperatorAQ::~ObsOperatorAQ() {}

// -----------------------------------------------------------------------------

void ObsOperatorAQ::simulateObs(const GomAQ & gvals, ObsVecAQ & yy, const ObsBias & bias,
                                ObsVecAQ &, ObsDiagsAQ &) const {
  oper_->simulateObs(gvals, yy, bias);
}

// -----------------------------------------------------------------------------

const oops::Variables & ObsOperatorAQ::requiredVars() const {
  return oper_->requiredVars();
}

// -----------------------------------------------------------------------------

std::unique_ptr<LocationsAQ> ObsOperatorAQ::locations() const {
  return oper_->locations();
}

// -----------------------------------------------------------------------------

void ObsOperatorAQ::print(std::ostream & os) const {
  os << *oper_;
}

// -----------------------------------------------------------------------------

}  // namespace aq
