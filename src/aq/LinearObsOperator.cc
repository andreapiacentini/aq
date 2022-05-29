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
