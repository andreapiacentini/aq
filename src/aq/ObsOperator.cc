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
