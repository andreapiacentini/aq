/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "aq/instantiateAqChangeVarFactory.h"
#include "aq/Traits.h"
#include "oops/runs/Run.h"
#include "oops/runs/Variational.h"
#include "saber/oops/instantiateCovarFactory.h"
#include "saber/oops/instantiateLocalizationFactory.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  aq::instantiateAqChangeVarFactory();
  saber::instantiateCovarFactory<aq::Traits>();
  saber::instantiateLocalizationFactory<aq::Traits>();
  oops::Variational<aq::Traits, aq::ObsTraits> var;
  return run.execute(var);
}
