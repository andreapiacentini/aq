/*
 * (C) Copyright 2017-2021 UCAR.
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_CHANGEVAR_H_
#define AQ_CHANGEVAR_H_

#include <ostream>
#include <string>

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
  static const std::string classname() {return "aq::ChangeVar";}

  ChangeVar(const eckit::Configuration &, const Geometry &);
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
