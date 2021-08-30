/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "aq/Future/ChangeVarTLADAQ.h"

#include <ostream>
#include <string>

#include "aq/GeometryAQ.h"
#include "aq/IncrementAQ.h"
#include "aq/StateAQ.h"
#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace aq {
// -----------------------------------------------------------------------------
ChangeVarTLADAQ::ChangeVarTLADAQ(const StateAQ &, const StateAQ &,
                                 const GeometryAQ & resol, const eckit::Configuration & conf) {}
// -----------------------------------------------------------------------------
ChangeVarTLADAQ::~ChangeVarTLADAQ() {}
// -----------------------------------------------------------------------------
void ChangeVarTLADAQ::multiply(const IncrementAQ & dxa, IncrementAQ & dxm) const {
  // AQ aq_change_var_tl_f90(dxa.fields().toFortran(), dxm.fields().toFortran());
  oops::Log::debug() << "ChangeVarTLADAQ::multiply" << dxm << std::endl;
}
// -----------------------------------------------------------------------------
void ChangeVarTLADAQ::multiplyInverse(const IncrementAQ & dxm, IncrementAQ & dxa) const {
  // AQ aq_change_var_tl_f90(dxm.fields().toFortran(), dxa.fields().toFortran());
  oops::Log::debug() << "ChangeVarTLADAQ::multiplyInverse" << dxm << std::endl;
}
// -----------------------------------------------------------------------------
void ChangeVarTLADAQ::multiplyAD(const IncrementAQ & dxm, IncrementAQ & dxa) const {
  // AQ aq_change_var_ad_f90(dxm.fields().toFortran(), dxa.fields().toFortran());
  oops::Log::debug() << "ChangeVarTLADAQ::multiplyAD" << dxm << std::endl;
}
// -----------------------------------------------------------------------------
void ChangeVarTLADAQ::multiplyInverseAD(const IncrementAQ & dxa, IncrementAQ & dxm) const {
  // AQ aq_change_var_ad_f90(dxa.fields().toFortran(), dxm.fields().toFortran());
  oops::Log::debug() << "ChangeVarTLADAQ::multiplyInverseAD" << dxm << std::endl;
}
// -----------------------------------------------------------------------------
void ChangeVarTLADAQ::print(std::ostream & os) const {
  os << "AQ linear change of variable";
}
// -----------------------------------------------------------------------------
}  // namespace aq

