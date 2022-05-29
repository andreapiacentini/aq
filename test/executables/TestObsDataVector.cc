/*
 * (C) Copyright 2021-2021 UCAR
 * (C) Copyright 2021-2022 CERFACS.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "aq/Traits.h"
#include "oops/runs/Run.h"
#include "test/interface/ObsDataVector.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::ObsDataVector<aq::ObsTraits> tests;
  return run.execute(tests);
}

