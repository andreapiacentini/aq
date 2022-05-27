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

#include "aq/ObsOperator.h"

#include "aq/GeoVals.h"
#include "aq/Locations.h"
#include "aq/ObsAuxControl.h"
#include "aq/ObsDiagnostics.h"
#include "aq/ObsOpBase.h"
#include "aq/ObsSpace.h"
#include "aq/ObsVec.h"
#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"

namespace aq {

// -----------------------------------------------------------------------------

ObsOperator::ObsOperator(const ObsSpace & os, const Parameters_ & params)
  : oper_(ObsOpFactory::create(os, params.config))
{}

// -----------------------------------------------------------------------------

ObsOperator::~ObsOperator() {}

// -----------------------------------------------------------------------------

void ObsOperator::simulateObs(const GeoVals & gvals, ObsVec & yy, const ObsAuxControl & bias,
                                ObsVec &, ObsDiagnostics &) const {
  oper_->simulateObs(gvals, yy, bias);
}

// -----------------------------------------------------------------------------

const oops::Variables & ObsOperator::requiredVars() const {
  return oper_->requiredVars();
}

// -----------------------------------------------------------------------------

std::unique_ptr<Locations> ObsOperator::locations() const {
  return oper_->locations();
}

// -----------------------------------------------------------------------------

void ObsOperator::print(std::ostream & os) const {
  os << *oper_;
}

// -----------------------------------------------------------------------------

}  // namespace aq
