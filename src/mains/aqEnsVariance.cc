/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
// AQ #include "aq/instantiateAqChangeVarFactory.h"
#include "aq/Traits.h"
#include "oops/runs/EnsVariance.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  // AQ aq::instantiateAqChangeVarFactory();
  oops::EnsVariance<aq::Traits> var;
  return run.execute(var);
}
