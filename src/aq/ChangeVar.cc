/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "aq/ChangeVar.h"

#include <ostream>
#include <string>

#include "oops/util/Logger.h"

#include "aq/AqFortran.h"
#include "aq/State.h"

namespace aq {
// -----------------------------------------------------------------------------
ChangeVar::ChangeVar(const Geometry &, const eckit::Configuration &) {}
// -----------------------------------------------------------------------------
ChangeVar::~ChangeVar() {}
// -----------------------------------------------------------------------------
void ChangeVar::changeVar(const State & xa, State & xm) const {
  xm = xa;
  // AQ  aq_change_var_f90(xa.fields().toFortran(), xm.fields().toFortran());
}
// -----------------------------------------------------------------------------
void ChangeVar::changeVarInverse(const State & xm, State & xa) const {
  xa = xm;
  // AQ aq_change_var_f90(xm.fields().toFortran(), xa.fields().toFortran());
}
// -----------------------------------------------------------------------------
void ChangeVar::print(std::ostream & os) const {
  os << "AQ change of variable";
}
// -----------------------------------------------------------------------------
}  // namespace aq


