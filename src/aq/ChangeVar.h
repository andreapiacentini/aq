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

#ifndef AQ_CHANGEVAR_H_
#define AQ_CHANGEVAR_H_

#include <ostream>
#include <string>

#include "aq/ChangeVarParameters.h"
#include "aq/Geometry.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace aq {
  class Geometry;
  class State;

// -----------------------------------------------------------------------------
/// AQ change of variable

class ChangeVar: public util::Printable {
 public:
  typedef ChangeVarParameters Parameters_;
  static const std::string classname() {return "aq::ChangeVar";}

  ChangeVar(const Parameters_ &, const Geometry &);
  ~ChangeVar();

/// Perform transforms
  void changeVar(State &, const oops::Variables &) const;
  void changeVarInverse(State &, const oops::Variables &) const;

 private:
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_CHANGEVAR_H_
