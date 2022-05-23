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

#ifndef AQ_LINEARCHANGEVAR_H_
#define AQ_LINEARCHANGEVAR_H_

#include <ostream>
#include <string>

#include "aq/LinearChangeVarParameters.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace aq {
  class Geometry;
  class State;
  class Increment;

// -----------------------------------------------------------------------------
/// AQ linear change of variable

class LinearChangeVar: public util::Printable {
 public:
  typedef LinearChangeVarParameters Parameters_;
  static const std::string classname() {return "aq::LinearChangeVar";}

  LinearChangeVar(const Geometry &, const Parameters_ &);
  ~LinearChangeVar();

/// Perform linear transforms
  void multiply(Increment &, const oops::Variables &) const;
  void multiplyInverse(Increment &, const oops::Variables &) const;
  void multiplyAD(Increment &, const oops::Variables &) const;
  void multiplyInverseAD(Increment &, const oops::Variables &) const;
  void setTrajectory(const State &, const State &);

 private:
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_LINEARCHANGEVAR_H_
