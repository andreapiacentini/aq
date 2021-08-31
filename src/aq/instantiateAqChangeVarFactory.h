/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_INSTANTIATEAQCHANGEVARFACTORY_H_
#define AQ_INSTANTIATEAQCHANGEVARFACTORY_H_

#include "oops/interface/LinearVariableChange.h"
#include "oops/interface/VariableChange.h"

#include "aq/AqTraits.h"
#include "aq/ChangeVarAQ.h"
#include "aq/ChangeVarTLADAQ.h"

namespace aq {

void instantiateAqChangeVarFactory() {
  static oops::VariableChangeMaker<AqTraits, ChangeVarAQ> makerchangevar_("ChVarAQ");
  static oops::VariableChangeMaker<AqTraits, ChangeVarAQ> makerdefchavar_("default");

  static oops::LinearVariableChangeMaker<aq::AqTraits,
                                   oops::LinearVariableChange<aq::AqTraits, aq::ChangeVarTLADAQ> >
               makerChLinVarAQ_("ChVarAQ");
}

}  // namespace aq

#endif  // AQ_INSTANTIATEAQCHANGEVARFACTORY_H_
