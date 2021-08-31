/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "aq/AqTraits.h"
#include "aq/instantiateAqChangeVarFactory.h"
#include "oops/runs/Run.h"
#include "saber/oops/instantiateCovarFactory.h"
#include "saber/oops/instantiateVariableChangeFactory.h"

#include "oops/runs/Variational.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  aq::instantiateAqChangeVarFactory();
  saber::instantiateCovarFactory<aq::AqTraits>();
  saber::instantiateLocalizationFactory<aq::AqTraits>();
  saber::instantiateVariableChangeFactory<aq::AqTraits>();
  oops::Variational<aq::AqTraits, aq::AqObsTraits> var;
  return run.execute(var);
}
