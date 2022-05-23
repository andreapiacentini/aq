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
