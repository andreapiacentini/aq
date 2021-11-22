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

#include "aq/ChangeVar.h"
#include "aq/LinearChangeVar.h"
#include "aq/Traits.h"

namespace aq {

void instantiateAqChangeVarFactory() {
  static oops::VariableChangeMaker<Traits, ChangeVar> makerchangevar_("ChVarAQ");
  static oops::VariableChangeMaker<Traits, ChangeVar> makerdefchavar_("default");

  static oops::LinearVariableChangeMaker<aq::Traits,
                                   oops::LinearVariableChange<aq::Traits, aq::LinearChangeVar> >
               makerChLinVarAQ_("ChVarAQ");
  static oops::LinearVariableChangeMaker<aq::Traits,
                                   oops::LinearVariableChange<aq::Traits, aq::LinearChangeVar> >
               makerChLinVarDef_("default");
}

}  // namespace aq

#endif  // AQ_INSTANTIATEAQCHANGEVARFACTORY_H_
