/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_FUTURE_CHANGEVARTLADAQ_H_
#define AQ_FUTURE_CHANGEVARTLADAQ_H_

#include <ostream>
#include <string>

#include "oops/util/Printable.h"

#include "oops/aq/AqFortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace aq {
  class GeometryAQ;
  class StateAQ;
  class IncrementAQ;

// -----------------------------------------------------------------------------
/// AQ linear change of variable

class ChangeVarTLADAQ: public util::Printable {
 public:
  static const std::string classname() {return "aq::ChangeVarAQ";}

  ChangeVarTLADAQ(const StateAQ &, const StateAQ &, const GeometryAQ &,
                  const eckit::Configuration &);
  ~ChangeVarTLADAQ();

/// Perform linear transforms
  void multiply(const IncrementAQ &, IncrementAQ &) const;
  void multiplyInverse(const IncrementAQ &, IncrementAQ &) const;
  void multiplyAD(const IncrementAQ &, IncrementAQ &) const;
  void multiplyInverseAD(const IncrementAQ &, IncrementAQ &) const;

 private:
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_FUTURE_CHANGEVARTLADAQ_H_
