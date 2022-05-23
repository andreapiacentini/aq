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
