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

#include "aq/LinearChangeVar.h"

#include <ostream>
#include <string>

#include "aq/Geometry.h"
#include "aq/Increment.h"
#include "aq/interface.h"
#include "aq/State.h"
#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace aq {
// -----------------------------------------------------------------------------
LinearChangeVar::LinearChangeVar(const Geometry &, const Parameters_ &) {}
// -----------------------------------------------------------------------------
LinearChangeVar::~LinearChangeVar() {}
// -----------------------------------------------------------------------------
void LinearChangeVar::multiply(Increment & dx, const oops::Variables & vars) const {
// AQ  aq_change_var_tl_f90(dx.fields().toFortran(), vars);
}
// -----------------------------------------------------------------------------
void LinearChangeVar::multiplyInverse(Increment & dx, const oops::Variables & vars) const {
// AQ  aq_change_var_tl_f90(dx.fields().toFortran(), vars);
}
// -----------------------------------------------------------------------------
void LinearChangeVar::multiplyAD(Increment & dx, const oops::Variables & vars) const {
// AQ  aq_change_var_ad_f90(dx.fields().toFortran(), vars);
}
// -----------------------------------------------------------------------------
void LinearChangeVar::multiplyInverseAD(Increment & dx, const oops::Variables & vars) const {
// AQ  aq_change_var_ad_f90(dx.fields().toFortran(), vars);
}
// -----------------------------------------------------------------------------
void LinearChangeVar::setTrajectory(const State & background, const State & firstGuess) {
  // No AQ trajectory used. No fortran to call here.
}
// -----------------------------------------------------------------------------
void LinearChangeVar::print(std::ostream & os) const {
  os << "AQ linear change of variable";
}
// -----------------------------------------------------------------------------
}  // namespace aq
