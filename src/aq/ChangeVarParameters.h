/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef AQ_CHANGEVARPARAMETERS_H_
#define AQ_CHANGEVARPARAMETERS_H_

#include <ostream>
#include <string>

#include "oops/base/VariableChangeParametersBase.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/Printable.h"

namespace aq {

// -------------------------------------------------------------------------------------------------
/// No change of variable parameters

class ChangeVarParameters : public oops::VariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(ChangeVarParameters, VariableChangeParametersBase)
 public:
};

// -------------------------------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_CHANGEVARPARAMETERS_H_
