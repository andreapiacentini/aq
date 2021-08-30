/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsOperatorTLAD.h"

#include "eckit/config/Configuration.h"
#include "model/GomAQ.h"
#include "model/ObsBias.h"
#include "model/ObsBiasIncrement.h"
#include "model/ObsOpBaseTLAD.h"
#include "model/ObsSpaceAQ.h"
#include "model/ObsVecAQ.h"
#include "oops/base/Variables.h"

namespace aq {

// -----------------------------------------------------------------------------

ObsOperatorTLAD::ObsOperatorTLAD(const ObsSpaceAQ & os, const eckit::Configuration & conf)
  : oper_(ObsOpTLADFactory::create(os, conf))
{}

// -----------------------------------------------------------------------------

ObsOperatorTLAD::~ObsOperatorTLAD() {}

// -----------------------------------------------------------------------------

void ObsOperatorTLAD::setTrajectory(const GomAQ & gvals, const ObsBias & bias) {
  oper_->setTrajectory(gvals, bias);
}

// -----------------------------------------------------------------------------

void ObsOperatorTLAD::simulateObsTL(const GomAQ & gvals, ObsVecAQ & yy,
                                    const ObsBiasIncrement & bias) const {
  oper_->simulateObsTL(gvals, yy, bias);
}

// -----------------------------------------------------------------------------

void ObsOperatorTLAD::simulateObsAD(GomAQ & gvals, const ObsVecAQ & yy,
                                    ObsBiasIncrement & bias) const {
  oper_->simulateObsAD(gvals, yy, bias);
}

// -----------------------------------------------------------------------------

const oops::Variables & ObsOperatorTLAD::requiredVars() const {
  return oper_->requiredVars();
}

// -----------------------------------------------------------------------------

void ObsOperatorTLAD::print(std::ostream & os) const {
  os << *oper_;
}

// -----------------------------------------------------------------------------

}  // namespace aq
