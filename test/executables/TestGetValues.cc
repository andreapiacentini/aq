/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/AqTraits.h"
#include "oops/runs/Run.h"
#include "test/interface/GetValues.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::GetValues<aq::AqTraits, aq::AqObsTraits> tests;
  return run.execute(tests);
}
