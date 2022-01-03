/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "aq/ChangeVar.h"

#include <ostream>
#include <string>

#include "oops/util/abor1_cpp.h"

#include "aq/State.h"

namespace aq {
// -----------------------------------------------------------------------------
ChangeVar::ChangeVar(const Parameters_ & params, const Geometry & geom) {}
// -----------------------------------------------------------------------------
ChangeVar::~ChangeVar() {}
// -----------------------------------------------------------------------------
void ChangeVar::changeVar(State & xx, const oops::Variables & vars) const {
  // AQ  aq_change_var_f90(xx.fields().toFortran(), vars);
}
// -----------------------------------------------------------------------------
void ChangeVar::changeVarInverse(State & xx, const oops::Variables & vars) const {
  // AQ  aq_change_var_f90(xx.fields().toFortran(), vars);
}
// -----------------------------------------------------------------------------
void ChangeVar::print(std::ostream & os) const {
  os << "AQ change of variable";
}
// -----------------------------------------------------------------------------
}  // namespace aq
