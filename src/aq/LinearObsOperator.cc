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

#include "aq/LinearObsOperator.h"

#include "aq/GeoVals.h"
#include "aq/LinearObsOpBase.h"
#include "aq/ObsAuxControl.h"
#include "aq/ObsAuxIncrement.h"
#include "aq/ObsSpace.h"
#include "aq/ObsVec.h"
#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"

namespace aq {

// -----------------------------------------------------------------------------

LinearObsOperator::LinearObsOperator(const ObsSpace & os, const Parameters_ & params)
  : oper_(LinearObsOpFactory::create(os, params.config))
{}

// -----------------------------------------------------------------------------

LinearObsOperator::~LinearObsOperator() {}

// -----------------------------------------------------------------------------

void LinearObsOperator::setTrajectory(const GeoVals & gvals, const ObsAuxControl & bias) {
  oper_->setTrajectory(gvals, bias);
}

// -----------------------------------------------------------------------------

void LinearObsOperator::simulateObsTL(const GeoVals & gvals, ObsVec & yy,
                                    const ObsAuxIncrement & bias) const {
  oper_->simulateObsTL(gvals, yy, bias);
}

// -----------------------------------------------------------------------------

void LinearObsOperator::simulateObsAD(GeoVals & gvals, const ObsVec & yy,
                                    ObsAuxIncrement & bias) const {
  oper_->simulateObsAD(gvals, yy, bias);
}

// -----------------------------------------------------------------------------

const oops::Variables & LinearObsOperator::requiredVars() const {
  return oper_->requiredVars();
}

// -----------------------------------------------------------------------------

void LinearObsOperator::print(std::ostream & os) const {
  os << *oper_;
}

// -----------------------------------------------------------------------------

}  // namespace aq
