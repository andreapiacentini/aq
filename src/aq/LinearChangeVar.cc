/*
 * (C) Copyright 2017-2018  UCAR.
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
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
LinearChangeVar::LinearChangeVar(const Geometry &, const eckit::Configuration &) {}
// -----------------------------------------------------------------------------
LinearChangeVar::~LinearChangeVar() {}
// -----------------------------------------------------------------------------
void LinearChangeVar::changeVarTL(Increment & dx, const oops::Variables & vars) const {
// AQ  aq_change_var_tl_f90(dx.fields().toFortran(), vars);
}
// -----------------------------------------------------------------------------
void LinearChangeVar::changeVarInverseTL(Increment & dx, const oops::Variables & vars) const {
// AQ  aq_change_var_tl_f90(dx.fields().toFortran(), vars);
}
// -----------------------------------------------------------------------------
void LinearChangeVar::changeVarAD(Increment & dx, const oops::Variables & vars) const {
// AQ  aq_change_var_ad_f90(dx.fields().toFortran(), vars);
}
// -----------------------------------------------------------------------------
void LinearChangeVar::changeVarInverseAD(Increment & dx, const oops::Variables & vars) const {
// AQ  aq_change_var_ad_f90(dx.fields().toFortran(), vars);
}
// -----------------------------------------------------------------------------
  void LinearChangeVar::changeVarTraj(const State & firstGuess, const oops::Variables & vars) {
  // No AQ trajectory used. No fortran to call here.
}
// -----------------------------------------------------------------------------
void LinearChangeVar::print(std::ostream & os) const {
  os << "AQ linear change of variable";
}
// -----------------------------------------------------------------------------
}  // namespace aq
