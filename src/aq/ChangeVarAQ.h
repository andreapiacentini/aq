/*
 * (C) Copyright 2017-2021 UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_CHANGEVARAQ_H_
#define AQ_CHANGEVARAQ_H_

#include <ostream>
#include <string>

#include "oops/base/VariableChangeBase.h"

#include "aq/AqTraits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace aq {
  class GeometryAQ;
  class StateAQ;

// -----------------------------------------------------------------------------
/// AQ change of variable

class ChangeVarAQ: public oops::VariableChangeBase<AqTraits> {
 public:
  static const std::string classname() {return "aq::ChangeVarAQ";}

  ChangeVarAQ(const GeometryAQ &, const eckit::Configuration &);
  ~ChangeVarAQ();

/// Perform transforms
  void changeVar(const StateAQ &, StateAQ &) const override;
  void changeVarInverse(const StateAQ &, StateAQ &) const override;

 private:
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_CHANGEVARAQ_H_
