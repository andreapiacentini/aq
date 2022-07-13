/*
 * (C) Copyright 2017-2018  UCAR.
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
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
  void changeVarTL(Increment &, const oops::Variables &) const;
  void changeVarInverseTL(Increment &, const oops::Variables &) const;
  void changeVarAD(Increment &, const oops::Variables &) const;
  void changeVarInverseAD(Increment &, const oops::Variables &) const;
  void changeVarTraj(const State &, const oops::Variables &);

 private:
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_LINEARCHANGEVAR_H_
