/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_LINEARCHANGEVARAQ_H_
#define AQ_LINEARCHANGEVARAQ_H_

#include <ostream>
#include <string>

#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace aq {
  class Geometry;
  class State;
  class Increment;

// -----------------------------------------------------------------------------
/// AQ linear change of variable

class LinearChangeVar: public util::Printable {
 public:
  static const std::string classname() {return "aq::ChangeVar";}

  LinearChangeVar(const State &, const State &, const Geometry &,
                  const eckit::Configuration &);
  ~LinearChangeVar();

/// Perform linear transforms
  void multiply(const Increment &, Increment &) const;
  void multiplyInverse(const Increment &, Increment &) const;
  void multiplyAD(const Increment &, Increment &) const;
  void multiplyInverseAD(const Increment &, Increment &) const;

 private:
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_LINEARCHANGEVARAQ_H_
