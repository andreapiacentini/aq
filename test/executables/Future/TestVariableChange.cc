/*
 * (C) Copyright 2020-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/Traits.h"
#include "model/instantiateAqChangeVarFactory.h"
#include "oops/runs/Run.h"
#include "test/interface/VariableChange.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  aq::instantiateAqChangeVarFactory();
  test::VariableChange<aq::Traits> tests;
  return run.execute(tests);
}