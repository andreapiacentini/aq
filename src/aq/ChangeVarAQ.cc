/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "aq/ChangeVarAQ.h"

#include <ostream>
#include <string>

#include "oops/util/Logger.h"

#include "aq/AqFortran.h"
#include "aq/StateAQ.h"

namespace aq {
// -----------------------------------------------------------------------------
ChangeVarAQ::ChangeVarAQ(const GeometryAQ &, const eckit::Configuration &) {}
// -----------------------------------------------------------------------------
ChangeVarAQ::~ChangeVarAQ() {}
// -----------------------------------------------------------------------------
void ChangeVarAQ::changeVar(const StateAQ & xa, StateAQ & xm) const {
  xm = xa;
  // AQ  aq_change_var_f90(xa.fields().toFortran(), xm.fields().toFortran());
}
// -----------------------------------------------------------------------------
void ChangeVarAQ::changeVarInverse(const StateAQ & xm, StateAQ & xa) const {
  xa = xm;
  // AQ aq_change_var_f90(xm.fields().toFortran(), xa.fields().toFortran());
}
// -----------------------------------------------------------------------------
void ChangeVarAQ::print(std::ostream & os) const {
  os << "AQ change of variable";
}
// -----------------------------------------------------------------------------
}  // namespace aq


