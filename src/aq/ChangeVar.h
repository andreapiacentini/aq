/*
 * (C) Copyright 2017-2021 UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_CHANGEVAR_H_
#define AQ_CHANGEVAR_H_

#include <ostream>
#include <string>

#include "oops/base/VariableChangeBase.h"

#include "aq/Traits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace aq {
  class Geometry;
  class State;

// -----------------------------------------------------------------------------
/// AQ change of variable

class ChangeVar: public oops::VariableChangeBase<Traits> {
 public:
  static const std::string classname() {return "aq::ChangeVar";}

  ChangeVar(const Geometry &, const eckit::Configuration &);
  ~ChangeVar();

/// Perform transforms
  void changeVar(const State &, State &) const override;
  void changeVarInverse(const State &, State &) const override;

 private:
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_CHANGEVAR_H_
