/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "aq/LinearChangeVar.h"

#include <ostream>
#include <string>

#include "aq/Geometry.h"
#include "aq/Increment.h"
#include "aq/State.h"
#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace aq {
// -----------------------------------------------------------------------------
LinearChangeVar::LinearChangeVar(const State &, const State &,
                                 const Geometry & resol, const eckit::Configuration & conf) {}
// -----------------------------------------------------------------------------
LinearChangeVar::~LinearChangeVar() {}
// -----------------------------------------------------------------------------
void LinearChangeVar::multiply(const Increment & dxa, Increment & dxm) const {
  dxm = dxa;
  // AQ aq_change_var_tl_f90(dxa.fields().toFortran(), dxm.fields().toFortran());
  oops::Log::debug() << "LinearChangeVar::multiply" << dxm << std::endl;
}
// -----------------------------------------------------------------------------
void LinearChangeVar::multiplyInverse(const Increment & dxm, Increment & dxa) const {
  dxa = dxm;
  // AQ aq_change_var_tl_f90(dxm.fields().toFortran(), dxa.fields().toFortran());
  oops::Log::debug() << "LinearChangeVar::multiplyInverse" << dxm << std::endl;
}
// -----------------------------------------------------------------------------
void LinearChangeVar::multiplyAD(const Increment & dxm, Increment & dxa) const {
  dxa = dxm;
  // AQ aq_change_var_ad_f90(dxm.fields().toFortran(), dxa.fields().toFortran());
  oops::Log::debug() << "LinearChangeVar::multiplyAD" << dxm << std::endl;
}
// -----------------------------------------------------------------------------
void LinearChangeVar::multiplyInverseAD(const Increment & dxa, Increment & dxm) const {
  dxm = dxa;
  // AQ aq_change_var_ad_f90(dxa.fields().toFortran(), dxm.fields().toFortran());
  oops::Log::debug() << "LinearChangeVar::multiplyInverseAD" << dxm << std::endl;
}
// -----------------------------------------------------------------------------
void LinearChangeVar::print(std::ostream & os) const {
  os << "AQ linear change of variable";
}
// -----------------------------------------------------------------------------
}  // namespace aq

